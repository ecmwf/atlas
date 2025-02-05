/*
 * (C) Copyright 2013 ECMWF
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <string>
#include <vector>

#include "atlas/functionspace/PointCloud.h"
#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/option/Options.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Point.h"
#include "atlas/util/Unique.h"


#include "eckit/mpi/Comm.h"
#include "eckit/log/Bytes.h"

#if ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif

namespace atlas {
namespace functionspace {

namespace detail {

namespace {

template <typename T, typename Field>
array::LocalView<T, 3> make_leveled_view(Field& field) {
    using namespace array;
    if (field.levels()) {
        if (field.variables()) {
            return make_view<T, 3>(field).slice(Range::all(), Range::all(), Range::all());
        }
        else {
            return make_view<T, 2>(field).slice(Range::all(), Range::all(), Range::dummy());
        }
    }
    else {
        if (field.variables()) {
            return make_view<T, 2>(field).slice(Range::all(), Range::dummy(), Range::all());
        }
        else {
            return make_view<T, 1>(field).slice(Range::all(), Range::dummy(), Range::dummy());
        }
    }
}

}  // namespace


static std::string get_mpi_comm(const eckit::Configuration& config) {
    if(config.has("mpi_comm")) {
        return config.getString("mpi_comm");
    }
    return mpi::comm().name();
}

template <>
PointCloud::PointCloud(const std::vector<PointXY>& points, const eckit::Configuration& config) {
    mpi_comm_ = get_mpi_comm(config);
    lonlat_     = Field("lonlat", array::make_datatype<double>(), array::make_shape(points.size(), 2));
    auto lonlat = array::make_view<double, 2>(lonlat_);
    for (idx_t j = 0, size = points.size(); j < size; ++j) {
        lonlat(j, 0) = points[j].x();
        lonlat(j, 1) = points[j].y();
    }
}

template <>
PointCloud::PointCloud(const std::vector<PointXYZ>& points, const eckit::Configuration& config) {
    mpi_comm_ = get_mpi_comm(config);
    lonlat_       = Field("lonlat", array::make_datatype<double>(), array::make_shape(points.size(), 2));
    vertical_     = Field("vertical", array::make_datatype<double>(), array::make_shape(points.size()));
    auto lonlat   = array::make_view<double, 2>(lonlat_);
    auto vertical = array::make_view<double, 1>(vertical_);
    for (idx_t j = 0, size = points.size(); j < size; ++j) {
        lonlat(j, 0) = points[j].x();
        lonlat(j, 1) = points[j].y();
        vertical(j)  = points[j].z();
    }
}

PointCloud::PointCloud(const Field& lonlat, const eckit::Configuration& config): lonlat_(lonlat) {
        mpi_comm_ = get_mpi_comm(config);
}

PointCloud::PointCloud(const Field& lonlat, const Field& ghost, const eckit::Configuration& config): lonlat_(lonlat), ghost_(ghost) {
    mpi_comm_ = get_mpi_comm(config);
    setupHaloExchange();
    setupGatherScatter();
}

PointCloud::PointCloud(const FieldSet& flds, const eckit::Configuration& config): lonlat_(flds["lonlat"]) {
    mpi_comm_ = get_mpi_comm(config);
    if (flds.has("ghost")) {
        ghost_ = flds["ghost"];
    }
    if (flds.has("remote_index")) {
        remote_index_ = flds["remote_index"];
    }
    if (flds.has("partition")) {
        partition_ = flds["partition"];
    }
    if (flds.has("global_index")) {
        global_index_ = flds["global_index"];
    }
    if( ghost_ && remote_index_ && partition_ ) {
        setupHaloExchange();
        setupGatherScatter();
    }
}

grid::Partitioner make_partitioner(const Grid& grid, const eckit::Configuration& config) {
    auto mpi_comm = get_mpi_comm(config);
    auto partitioner = grid.partitioner();
    if( config.has("partitioner") ) {
        partitioner.set("type",config.getString("partitioner"));
    }
    if( not partitioner.has("type") ) {
        partitioner.set("type","equal_regions");
    }
    partitioner.set("mpi_comm",mpi_comm);
    return grid::Partitioner(partitioner);
}

PointCloud::PointCloud(const Grid& grid, const eckit::Configuration& config) :
    PointCloud(grid, make_partitioner(grid,config), config) {
}

PointCloud::PointCloud(const Grid& grid, const grid::Partitioner& _partitioner, const eckit::Configuration& config) {
    ATLAS_TRACE("PointCloud(grid,partitioner,config)");
    mpi_comm_ = get_mpi_comm(config);
    auto& comm = mpi::comm(mpi_comm_);
    double halo_radius;
    config.get("halo_radius", halo_radius = 0.);

    grid::Partitioner partitioner(_partitioner);
    if ( not partitioner ) {
        partitioner = grid::Partitioner("equal_regions", util::Config("mpi_comm",mpi_comm_));
    }
    part_ = comm.rank();

    nb_partitions_ = partitioner.nb_partitions();
    auto distribution = partitioner.partition(grid);
    auto size_owned = distribution.nb_pts()[part_];
    size_owned_ = size_owned;

    if (halo_radius == 0. || nb_partitions_ == 1) {
        idx_t size_halo = size_owned;
        ATLAS_ASSERT(size_owned > 0);
        lonlat_       = Field("lonlat", array::make_datatype<double>(), array::make_shape(size_halo, 2));
        ghost_        = Field("ghost", array::make_datatype<int>(), array::make_shape(size_halo));
        global_index_ = Field("global_index", array::make_datatype<gidx_t>(), array::make_shape(size_halo));
        partition_    = Field("partition", array::make_datatype<int>(), array::make_shape(size_halo));
        remote_index_ = Field("remote_index", array::make_datatype<idx_t>(), array::make_shape(size_halo));

        auto lonlat = array::make_view<double, 2>(lonlat_);
        auto ridx = array::make_indexview<idx_t,1>(remote_index_);
        auto gidx = array::make_view<gidx_t,1>(global_index_);
        array::make_view<int,1>(ghost_).assign(0);
        array::make_view<int,1>(partition_).assign(part_);

        idx_t j{0};
        gidx_t g{0};
        for (auto p : grid.lonlat()) {
            if( distribution.partition(g) == part_ ) {
                gidx(j) = g+1;
                ridx(j) = j;
                lonlat(j, 0) = p.lon();
                lonlat(j, 1) = p.lat();
                ++j;
            }
            ++g;
        }
    }
    else {

        std::vector<int> keep;

        {
            ATLAS_TRACE("Build list of points to keep");
            std::vector<PointLonLat> owned_lonlat;
            owned_lonlat.reserve(size_owned);

            auto kdtree = util::IndexKDTree(config);
            {
                ATLAS_TRACE("build kdtree");
                kdtree.reserve(grid.size());
                idx_t j{0};
                for (auto p : grid.lonlat()) {
                    if( distribution.partition(j) == part_ ) {
                        owned_lonlat.emplace_back(p);
                    }
                    kdtree.insert(p,j);
                    ++j;
                }
                kdtree.build();
            }

            keep.resize(grid.size());
            {
                ATLAS_TRACE("search kdtree");
                for (idx_t j=0; j<size_owned; ++j) {
                    const auto& p = owned_lonlat[j];
                    auto points = kdtree.closestPointsWithinRadius(p,halo_radius).payloads();
                    for( idx_t jj: points ) {
                        keep[jj] = 1;
                    }
                }
            }
        }

        {
            ATLAS_TRACE("create fields");

            auto size_halo = std::accumulate(keep.begin(),keep.end(),0);

            lonlat_         = Field("lonlat", array::make_datatype<double>(), array::make_shape(size_halo, 2));
            partition_      = Field("partition", array::make_datatype<int>(), array::make_shape(size_halo));
            ghost_          = Field("ghost", array::make_datatype<int>(), array::make_shape(size_halo));
            global_index_   = Field("global_index", array::make_datatype<gidx_t>(), array::make_shape(size_halo));
            max_glb_idx_    = grid.size();
            auto lonlat     = array::make_view<double,2>(lonlat_);
            auto partition  = array::make_view<int,1>(partition_);
            auto ghost      = array::make_view<int,1>(ghost_);
            auto glb_idx    = array::make_view<gidx_t,1>(global_index_);

            gidx_t g=0;
            idx_t j = 0;
            for (auto p : grid.lonlat()) {
                if (keep[g]) {
                    lonlat(j, 0) = p.lon();
                    lonlat(j, 1) = p.lat();
                    partition(j) = distribution.partition(g);
                    ghost(j) = partition(j) != part_;
                    glb_idx(j) = g+1;
                    ++j;
                }
                ++g;
            }

        }

    }

    setupHaloExchange();
    setupGatherScatter();
}


Field PointCloud::ghost() const {
    if (not ghost_) {
        ghost_ = Field("ghost", array::make_datatype<int>(), array::make_shape(size()));
        array::make_view<int, 1>(ghost_).assign(0);
    }
    return ghost_;
}

array::ArrayShape PointCloud::config_shape(const eckit::Configuration& config) const {
    idx_t _size  = size();
    bool global(false);
    if (config.get("global", global)) {
        if (global) {
            idx_t owner(0);
            config.get("owner", owner);
            idx_t rank = mpi::comm(mpi_comm()).rank();
            _size = (rank == owner ? size_global_ : 0);
        }
    }

    array::ArrayShape shape;

    shape.emplace_back(_size);

    idx_t levels(levels_);
    config.get("levels", levels);
    if (levels > 0) {
        shape.emplace_back(levels);
    }

    idx_t variables(0);
    config.get("variables", variables);
    if (variables > 0) {
        shape.emplace_back(variables);
    }

    return shape;
}

array::ArrayAlignment PointCloud::config_alignment(const eckit::Configuration& config) const {
    int alignment(1);
    config.get("alignment", alignment);
    return alignment;
}

array::ArraySpec PointCloud::config_spec(const eckit::Configuration& config) const {
    return array::ArraySpec(config_shape(config), config_alignment(config));
}

array::DataType PointCloud::config_datatype(const eckit::Configuration& config) const {
    array::DataType::kind_t kind;
    if (!config.get("datatype", kind)) {
        throw_Exception("datatype missing", Here());
    }
    return array::DataType(kind);
}

std::string PointCloud::config_name(const eckit::Configuration& config) const {
    std::string name;
    config.get("name", name);
    return name;
}

const parallel::HaloExchange& PointCloud::halo_exchange() const {
    return *halo_exchange_;
}

void PointCloud::gather(const FieldSet& local_fieldset, FieldSet& global_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());

    for (idx_t f = 0; f < local_fieldset.size(); ++f) {
        const Field& loc      = local_fieldset[f];
        Field& glb            = global_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root(0);
        glb.metadata().get("owner", root);

        if (loc.datatype() == array::DataType::kind<int>()) {
            parallel::Field<int const> loc_field(make_leveled_view<const int>(loc));
            parallel::Field<int> glb_field(make_leveled_view<int>(glb));
            gather().gather(&loc_field, &glb_field, nb_fields, root);
        }
        else if (loc.datatype() == array::DataType::kind<long>()) {
            parallel::Field<long const> loc_field(make_leveled_view<const long>(loc));
            parallel::Field<long> glb_field(make_leveled_view<long>(glb));
            gather().gather(&loc_field, &glb_field, nb_fields, root);
        }
        else if (loc.datatype() == array::DataType::kind<float>()) {
            parallel::Field<float const> loc_field(make_leveled_view<const float>(loc));
            parallel::Field<float> glb_field(make_leveled_view<float>(glb));
            gather().gather(&loc_field, &glb_field, nb_fields, root);
        }
        else if (loc.datatype() == array::DataType::kind<double>()) {
            parallel::Field<double const> loc_field(make_leveled_view<const double>(loc));
            parallel::Field<double> glb_field(make_leveled_view<double>(glb));
            gather().gather(&loc_field, &glb_field, nb_fields, root);
        }
        else {
            throw_Exception("datatype not supported", Here());
        }
    }
}

void PointCloud::gather(const Field& local, Field& global) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields, global_fields);
}
const parallel::GatherScatter& PointCloud::gather() const {
    ATLAS_ASSERT(gather_scatter_);
    return *gather_scatter_;
}
const parallel::GatherScatter& PointCloud::scatter() const {
    ATLAS_ASSERT(gather_scatter_);
    return *gather_scatter_;
}

void PointCloud::scatter(const FieldSet& global_fieldset, FieldSet& local_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());

    for (idx_t f = 0; f < local_fieldset.size(); ++f) {
        const Field& glb      = global_fieldset[f];
        Field& loc            = local_fieldset[f];
        const idx_t nb_fields = 1;
        idx_t root(0);
        glb.metadata().get("owner", root);

        if (loc.datatype() == array::DataType::kind<int>()) {
            parallel::Field<int const> glb_field(make_leveled_view<const int>(glb));
            parallel::Field<int> loc_field(make_leveled_view<int>(loc));
            scatter().scatter(&glb_field, &loc_field, nb_fields, root);
        }
        else if (loc.datatype() == array::DataType::kind<long>()) {
            parallel::Field<long const> glb_field(make_leveled_view<const long>(glb));
            parallel::Field<long> loc_field(make_leveled_view<long>(loc));
            scatter().scatter(&glb_field, &loc_field, nb_fields, root);
        }
        else if (loc.datatype() == array::DataType::kind<float>()) {
            parallel::Field<float const> glb_field(make_leveled_view<const float>(glb));
            parallel::Field<float> loc_field(make_leveled_view<float>(loc));
            scatter().scatter(&glb_field, &loc_field, nb_fields, root);
        }
        else if (loc.datatype() == array::DataType::kind<double>()) {
            parallel::Field<double const> glb_field(make_leveled_view<const double>(glb));
            parallel::Field<double> loc_field(make_leveled_view<double>(loc));
            scatter().scatter(&glb_field, &loc_field, nb_fields, root);
        }
        else {
            throw_Exception("datatype not supported", Here());
        }

        auto name = loc.name();
        glb.metadata().broadcast(loc.metadata(), root);
        loc.metadata().set("global", false);
        if( !name.empty() ) {
            loc.metadata().set("name", name);
        }
    }
}

void PointCloud::scatter(const Field& global, Field& local) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}


void PointCloud::set_field_metadata(const eckit::Configuration& config, Field& field) const {
    field.set_functionspace(this);

    bool global(false);
    if (config.get("global", global)) {
        if (global) {
            idx_t owner(0);
            config.get("owner", owner);
            field.metadata().set("owner", owner);
        }
    }
    field.metadata().set("global", global);

    idx_t levels(levels_);
    config.get("levels", levels);
    field.set_levels(levels);

    idx_t variables(0);
    config.get("variables", variables);
    field.set_variables(variables);

    if (config.has("type")) {
        field.metadata().set("type", config.getString("type"));
    }
}

Field PointCloud::createField(const eckit::Configuration& options) const {
    Field field(config_name(options), config_datatype(options), config_spec(options));
    set_field_metadata(options, field);
    return field;
}

Field PointCloud::createField(const Field& other, const eckit::Configuration& config) const {
    return createField(option::datatype(other.datatype()) | option::levels(other.levels()) |
                       option::variables(other.variables()) |
                       option::type(other.metadata().getString("type", "scalar")) | config);
}

std::string PointCloud::distribution() const {
    return (partition_ ?  std::string("pointcloud") : std::string("serial"));
}

static const array::Array& get_dummy() {
    static array::ArrayT<double> dummy_{1};
    return dummy_;
}

template <typename Point>
PointCloud::IteratorT<Point>::IteratorT(const atlas::functionspace::detail::PointCloud& fs, bool begin):
    fs_(fs),
    xy_(array::make_view<const double, 2>(fs_.lonlat())),
    z_(array::make_view<const double, 1>(bool(fs_.vertical()) ? fs_.vertical().array() : get_dummy())),
    n_(begin ? 0 : fs_.size()),
    size_(fs_.size()) {}

template <typename Point>
bool PointCloud::IteratorT<Point>::next(Point& p) {
    if (n_ < size_) {
        p[XX] = xy_(n_, XX);
        p[YY] = xy_(n_, YY);
        if (Point::DIMS == 3) {
            p[ZZ] = z_(n_);
        }
        ++n_;
        return true;
    }
    return false;
}

template <typename Point>
const Point PointCloud::IteratorT<Point>::operator*() const {
    Point p;
    p[XX] = xy_(n_, XX);
    p[YY] = xy_(n_, YY);
    if (Point::DIMS == 3) {
        p[ZZ] = z_(n_);
    }
    return p;
}

template class PointCloud::IteratorT<PointXYZ>;
template class PointCloud::IteratorT<PointXY>;
template class PointCloud::IteratorT<PointLonLat>;


namespace {

template <int RANK>
void dispatch_haloExchange(Field& field, const parallel::HaloExchange& halo_exchange, bool on_device) {
    if (field.datatype() == array::DataType::kind<int>()) {
        halo_exchange.template execute<int, RANK>(field.array(), on_device);
    }
    else if (field.datatype() == array::DataType::kind<long>()) {
        halo_exchange.template execute<long, RANK>(field.array(), on_device);
    }
    else if (field.datatype() == array::DataType::kind<float>()) {
        halo_exchange.template execute<float, RANK>(field.array(), on_device);
    }
    else if (field.datatype() == array::DataType::kind<double>()) {
        halo_exchange.template execute<double, RANK>(field.array(), on_device);
    }
    else {
        throw_Exception("datatype not supported", Here());
    }
    field.set_dirty(false);
}

template <int RANK>
void dispatch_adjointHaloExchange(Field& field, const parallel::HaloExchange& halo_exchange, bool on_device) {
    if (field.datatype() == array::DataType::kind<int>()) {
        halo_exchange.template execute_adjoint<int, RANK>(field.array(), on_device);
    }
    else if (field.datatype() == array::DataType::kind<long>()) {
        halo_exchange.template execute_adjoint<long, RANK>(field.array(), on_device);
    }
    else if (field.datatype() == array::DataType::kind<float>()) {
        halo_exchange.template execute_adjoint<float, RANK>(field.array(), on_device);
    }
    else if (field.datatype() == array::DataType::kind<double>()) {
        halo_exchange.template execute_adjoint<double, RANK>(field.array(), on_device);
    }
    else {
        throw_Exception("datatype not supported", Here());
    }
    field.set_dirty(false);
}

}  // namespace

void PointCloud::haloExchange(const FieldSet& fieldset, bool on_device) const {
    if (halo_exchange_) {
        for (idx_t f = 0; f < fieldset.size(); ++f) {
            Field& field = const_cast<FieldSet&>(fieldset)[f];
            switch (field.rank()) {
                case 1:
                    dispatch_haloExchange<1>(field, halo_exchange(), on_device);
                    break;
                case 2:
                    dispatch_haloExchange<2>(field, halo_exchange(), on_device);
                   break;
                case 3:
                    dispatch_haloExchange<3>(field, halo_exchange(), on_device);
                    break;
                case 4:
                    dispatch_haloExchange<4>(field, halo_exchange(), on_device);
                    break;
                default:
                    throw_Exception("Rank not supported", Here());
            }
            field.set_dirty(false);
        }
    }
}

void PointCloud::haloExchange(const Field& field, bool on_device) const {
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset, on_device);
}

void PointCloud::create_remote_index() const {
    ATLAS_TRACE();
    const auto& comm = mpi::comm(mpi_comm_);
    const int mpi_rank = comm.rank();
    const int mpi_size = comm.size();
    auto size_halo = lonlat_.shape(0);

    remote_index_ = Field("remote_idx", array::make_datatype<idx_t>(), array::make_shape(size_halo));
    auto remote_idx     = array::make_indexview<idx_t, 1>(remote_index_);
    auto ghost          = array::make_view<int, 1>(ghost_);

    atlas_omp_parallel_for(idx_t n = 0; n < size_halo; ++n) {
        if(not ghost(n)) {
            remote_idx(n) = n;
        }
        else {
            remote_idx(n) = -1;
        }
    }

    auto build_partition_graph = [this,mpi_size,mpi_rank,size_halo,&comm]() -> std::unique_ptr<Mesh::PartitionGraph> {

        auto part  = array::make_view<int, 1>(this->partition_);
        auto ghost = array::make_view<int, 1>(this->ghost_);

        std::vector<int> others_set(mpi_size, 0);
        others_set[mpi_rank] = 1;
        for( idx_t i = 0; i<size_halo; ++i) {
            if (ghost(i)) {
                others_set[part(i)] = 1;  // present
            }
        }
        std::vector<int> others;
        others.reserve(mpi_size);
        for (idx_t p = 0; p < mpi_size; ++p) {
            if (others_set[p]) {
                others.emplace_back(p);
            }
        }

        eckit::mpi::Buffer<int> recv_others(mpi_size);

        ATLAS_TRACE_MPI(ALLGATHER) { comm.allGatherv(others.begin(), others.end(), recv_others); }

        std::vector<idx_t> counts(recv_others.counts.begin(), recv_others.counts.end());
        std::vector<idx_t> displs(recv_others.displs.begin(), recv_others.displs.end());
        std::vector<idx_t> values(recv_others.buffer.begin(), recv_others.buffer.end());
        return std::unique_ptr<Mesh::PartitionGraph>(
            new Mesh::PartitionGraph(values.data(), mpi_size, displs.data(), counts.data()));
    };

    std::unique_ptr<Mesh::PartitionGraph> graph_ptr;
    ATLAS_TRACE_SCOPE("Building partition graph...") { graph_ptr = build_partition_graph(); }
    const Mesh::PartitionGraph& graph = *graph_ptr;

    ATLAS_TRACE_SCOPE("Setup remote_index fields...") {
        auto p = array::make_view<int, 1>(partition_);
        auto g = array::make_view<gidx_t, 1>(global_index_);

        auto neighbours           = graph.nearestNeighbours(mpi_rank);
        const idx_t nb_neighbours = static_cast<idx_t>(neighbours.size());
        std::unordered_map<int, idx_t> part_to_neighbour;
        ATLAS_TRACE_SCOPE("part_to_neighbour") {
            part_to_neighbour.reserve(nb_neighbours);
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                part_to_neighbour[neighbours[j]] = j;
            }
        }
        std::vector<idx_t> halo_per_neighbour(neighbours.size(), 0);
        ATLAS_TRACE_SCOPE("set halo_per_neighbour")
        for( idx_t i = 0; i<size_halo; ++i) {
            if (ghost(i)) {
                halo_per_neighbour[part_to_neighbour[p(i)]]++;
            }
        }


        std::vector<std::vector<gidx_t>> g_per_neighbour(neighbours.size());
        ATLAS_TRACE_SCOPE("assemble g_per_neighbour") {
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                g_per_neighbour[j].reserve(halo_per_neighbour[j]);
            }
            for( idx_t j = 0; j<size_halo; ++j) {
                if (ghost(j)) {
                    g_per_neighbour[part_to_neighbour[p(j)]].emplace_back(g(j));
                }
            }
        }

        std::vector<eckit::mpi::Request> send_requests(neighbours.size());
        std::vector<eckit::mpi::Request> recv_requests(neighbours.size());

        std::vector<idx_t> recv_size(neighbours.size());
        std::vector<idx_t> send_size(neighbours.size());

        int tag = 0;
        ATLAS_TRACE_SCOPE("send-receive g_per_neighbour size") {
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                send_size[j]     = static_cast<idx_t>(g_per_neighbour[j].size());
                send_requests[j] = comm.iSend(send_size[j], neighbours[j], tag);
                recv_requests[j] = comm.iReceive(recv_size[j], neighbours[j], tag);
            }

            ATLAS_TRACE_MPI(WAIT) {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    comm.wait(send_requests[j]);
                }

                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    comm.wait(recv_requests[j]);
                }
            }
        }

        std::vector<std::vector<gidx_t>> recv_g_per_neighbour(neighbours.size());
        ATLAS_TRACE_SCOPE("send-receive g_per_neighbour")
        for (idx_t j = 0; j < nb_neighbours; ++j) {
            recv_g_per_neighbour[j].resize(recv_size[j]);

            send_requests[j] = comm.iSend(g_per_neighbour[j].data(), g_per_neighbour[j].size(), neighbours[j], tag);
            recv_requests[j] =
                comm.iReceive(recv_g_per_neighbour[j].data(), recv_g_per_neighbour[j].size(), neighbours[j], tag);
        }

        std::vector<std::vector<idx_t>> send_r_per_neighbour(neighbours.size());

        // Assemble "send_r_per_neighbour"
        // Depending if we can afford memory for a globally sized array, we can have a
        // much faster version of g_to_r map using std::vector.
        // TODO: using c++14 we can create a polymorphic lambda to avoid duplicated
        // code in the two branches of following if.

        // size_t max_glb_idx = grid().size();
        ATLAS_ASSERT(max_glb_idx_);

        atlas::vector<idx_t> g_to_r_vector;
        bool use_unordered_map_fallback = false;
        try {
            g_to_r_vector.resize(max_glb_idx_ + 1);
        }
        catch (std::bad_alloc& e) {
            if (comm.size() > 1) {
                Log::warning() << "Could not allocate " << eckit::Bytes((max_glb_idx_ + 1) * sizeof(idx_t)) << Here()
                                << "\n"
                                << "Using slower unordered_map fallback to map global to remote indices"
                                << std::endl;
                use_unordered_map_fallback = true;
            }
            else {
                throw_Exception(
                    "Could not allocate " + std::string(eckit::Bytes((max_glb_idx_ + 1) * sizeof(idx_t))), Here());
            }
        }
        if (not use_unordered_map_fallback) {
            auto& g_to_r = g_to_r_vector;
            ATLAS_TRACE_SCOPE("g_to_r (using vector)") {
                atlas_omp_parallel_for(idx_t j = 0; j < size_halo; ++j) {
                    if (not ghost(j)) {
                        // parallel omp possible for ``` g_to_r[g(j)] = j ``` as we only loop over size_owned,
                        // where global_index is such that race-conditions cannot occur
#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
                        if (g(j) >= max_glb_idx_ + 1) {
                            ATLAS_DEBUG_VAR(g(j));
                            throw_OutOfRange("g_to_r", g(j), max_glb_idx_ + 1, Here());
                        }
#endif
                        g_to_r[g(j)] = j;
                    }
                }
            }
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                send_r_per_neighbour[j].reserve(recv_size[j]);

                comm.wait(recv_requests[j]);  // wait for recv_g_per_neighbour[j]
                for (gidx_t gidx : recv_g_per_neighbour[j]) {
                    send_r_per_neighbour[j].emplace_back(g_to_r[gidx]);
                }
            }
        }
        else {
            std::unordered_map<gidx_t, idx_t> g_to_r;
            g_to_r.reserve(size_halo);

            ATLAS_TRACE_SCOPE("g_to_r (using unordered set)") {
                for (idx_t j = 0; j < size_halo; ++j) {
                    if (not ghost(j)) {
                        g_to_r[g(j)] = j;
                    }
                }
            }
            ATLAS_TRACE_SCOPE("assemble r_per_neighbour")
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                send_r_per_neighbour[j].reserve(recv_size[j]);

                comm.wait(recv_requests[j]);  // wait for recv_g_per_neighbour[j]
                for (gidx_t gidx : recv_g_per_neighbour[j]) {
                    send_r_per_neighbour[j].emplace_back(g_to_r[gidx]);
                }
            }
        }

        std::vector<std::vector<idx_t>> r_per_neighbour(neighbours.size());

        ATLAS_TRACE_SCOPE("send-receive r_per_neighbour") {
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                r_per_neighbour[j].resize(halo_per_neighbour[j]);
            }

            ATLAS_TRACE_MPI(ALLTOALL) {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    comm.wait(send_requests[j]);
                    send_requests[j] = comm.iSend(send_r_per_neighbour[j].data(), send_r_per_neighbour[j].size(),
                                                    neighbours[j], tag);
                    recv_requests[j] =
                        comm.iReceive(r_per_neighbour[j].data(), r_per_neighbour[j].size(), neighbours[j], tag);
                }
            }
            ATLAS_TRACE_MPI(WAIT) {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    comm.wait(recv_requests[j]);
                }
            }
        }


        std::vector<idx_t> counters(neighbours.size(), 0);
        ATLAS_TRACE_SCOPE("set remote_idx")
        for (idx_t j = 0; j < size_halo; ++j) {
            if (ghost(j)) {
                idx_t neighbour = part_to_neighbour[p(j)];
                remote_idx(j)   = r_per_neighbour[neighbour][counters[neighbour]++];
            }
            ATLAS_ASSERT(remote_idx(j) >= 0);
        }

        ATLAS_TRACE_MPI(WAIT) {
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                comm.wait(send_requests[j]);
            }
        }
    }
}

void PointCloud::setupHaloExchange() {
    ATLAS_TRACE();
    if (ghost_ and partition_ and global_index_ and not remote_index_) {
        create_remote_index();
    }
    else if (not partition_ or not remote_index_) {
        ATLAS_TRACE("do setup");
        const auto& comm = mpi::comm(mpi_comm_);
        const int mpi_rank = comm.rank();
        const int mpi_size = comm.size();

        auto lonlat_v = array::make_view<double, 2>(lonlat_);
        // data structure containing a flag to identify the 'ghost points';
        // 0={is not a ghost point}, 1={is a ghost point}
        auto is_ghost = array::make_view<int, 1>(ghost_);

        std::vector<PointXY> opoints_local;
        std::vector<PointXY> gpoints_local;
        std::vector<uidx_t> lonlat_u;
        std::vector<uidx_t> opoints_local_u;
        lonlat_u.reserve(lonlat_v.shape(0));
        opoints_local_u.reserve(is_ghost.shape(0));
        gpoints_local.reserve(is_ghost.shape(0));
        opoints_local.reserve(is_ghost.shape(0));

        for (idx_t i = 0; i < lonlat_v.shape(0); ++i) {
            lonlat_u.emplace_back(util::unique_lonlat(lonlat_v(i, XX), lonlat_v(i, YY)));
        }

        idx_t j{0};
        for (idx_t i = 0; i < is_ghost.shape(0); ++i) {
            PointXY loc(lonlat_v(j, XX), lonlat_v(j, YY));
            if (is_ghost(i)) {
                gpoints_local.emplace_back(loc);
            }
            else {
                opoints_local.emplace_back(loc);
                opoints_local_u.emplace_back(util::unique_lonlat(loc.x(), loc.y()));
            }
            ++j;
        }

        std::vector<double> coords_gp_local;
        coords_gp_local.reserve(gpoints_local.size() * 2);

        for (auto& gp : gpoints_local) {
            coords_gp_local.push_back(gp[XX]);
            coords_gp_local.push_back(gp[YY]);
        }

        eckit::mpi::Buffer<double> buffers_rec(mpi_size);

        ATLAS_TRACE_SCOPE("mpi all gather") {
            comm.allGatherv(coords_gp_local.begin(), coords_gp_local.end(), buffers_rec);
        }

        std::vector<PointXY> gpoints_global;

        for (std::size_t pe = 0; pe < mpi_size; ++pe) {
            for (std::size_t j = 0; j < buffers_rec.counts[pe] / 2; ++j) {
                PointXY loc_gp(*(buffers_rec.begin() + buffers_rec.displs[pe] + 2 * j + XX),
                               *(buffers_rec.begin() + buffers_rec.displs[pe] + 2 * j + YY));
                gpoints_global.emplace_back(loc_gp);
            }
        }

        std::vector<uidx_t> gpoints_global_u;
        gpoints_global_u.reserve(gpoints_global.size());
        for (atlas::PointXY& loc : gpoints_global) {
            gpoints_global_u.emplace_back(util::unique_lonlat(loc.x(), loc.y()));
        }

        std::vector<int> partition_ids_gp_global(gpoints_global.size(), -1);
        std::vector<int> remote_index_gp_global(gpoints_global.size(), -1);

        std::vector<uidx_t>::iterator iter_xy_gp_01;

        ATLAS_TRACE_SCOPE("find global")
        for (std::size_t idx = 0; idx < gpoints_global_u.size(); ++idx) {
            iter_xy_gp_01 = std::find(opoints_local_u.begin(), opoints_local_u.end(), gpoints_global_u.at(idx));
            if (iter_xy_gp_01 != opoints_local_u.end()) {
                std::size_t ridx                = std::distance(opoints_local_u.begin(), iter_xy_gp_01);
                partition_ids_gp_global.at(idx) = mpi_rank;
                remote_index_gp_global.at(idx)  = ridx;
            }
        }

        ATLAS_TRACE_SCOPE("mpi all reduce") {
            comm.allReduceInPlace(partition_ids_gp_global.begin(), partition_ids_gp_global.end(), eckit::mpi::max());
            comm.allReduceInPlace(remote_index_gp_global.begin(), remote_index_gp_global.end(), eckit::mpi::max());
        }

        std::vector<int> partition_ids_local(lonlat_v.shape(0), -1);
        std::vector<idx_t> remote_index_local(lonlat_v.shape(0), -1);

        idx_t idx_loc{0};
        std::vector<uidx_t>::iterator iter_xy_gp_02;

        ATLAS_TRACE_SCOPE("find local")
        for (idx_t i = 0; i < lonlat_v.shape(0); ++i) {
            iter_xy_gp_02 = std::find(gpoints_global_u.begin(), gpoints_global_u.end(), lonlat_u.at(i));
            if (iter_xy_gp_02 != gpoints_global_u.end()) {
                std::size_t idx_gp           = std::distance(gpoints_global_u.begin(), iter_xy_gp_02);
                partition_ids_local[idx_loc] = partition_ids_gp_global[idx_gp];
                remote_index_local[idx_loc]  = remote_index_gp_global[idx_gp];
            }
            else {
                partition_ids_local[idx_loc] = mpi_rank;
                remote_index_local[idx_loc]  = idx_loc;
            }
            ++idx_loc;
        }

        partition_ = Field("partition", array::make_datatype<int>(), array::make_shape(partition_ids_local.size()));

        auto partitionv = array::make_view<int, 1>(partition_);
        for (idx_t i = 0; i < partitionv.shape(0); ++i) {
            partitionv(i) = partition_ids_local.at(i);
        }

        remote_index_ =
            Field("remote_index", array::make_datatype<idx_t>(), array::make_shape(remote_index_local.size()));

        auto remote_indexv = array::make_indexview<idx_t, 1>(remote_index_);
        for (idx_t i = 0; i < remote_indexv.shape(0); ++i) {
            remote_indexv(i) = remote_index_local.at(i);
        }
    }

    ATLAS_ASSERT(partition_);
    ATLAS_ASSERT(ghost_);
    ATLAS_ASSERT(remote_index_);
    ATLAS_ASSERT(ghost_.size() == remote_index_.size());
    ATLAS_ASSERT(ghost_.size() == partition_.size());
 
    halo_exchange_.reset(new parallel::HaloExchange());
    halo_exchange_->setup(mpi_comm_,
                          array::make_view<int, 1>(partition_).data(),
                          array::make_view<idx_t, 1>(remote_index_).data(),
                          REMOTE_IDX_BASE,
                          ghost_.size());
}

void PointCloud::setupGatherScatter() {
    if (ghost_ and partition_ and remote_index_ and global_index_) {
        gather_scatter_.reset(new parallel::GatherScatter());
        gather_scatter_->setup(mpi_comm_,
                            array::make_view<int, 1>(partition_).data(),
                            array::make_view<idx_t, 1>(remote_index_).data(),
                            REMOTE_IDX_BASE,
                            array::make_view<gidx_t, 1>(global_index_).data(),
                            array::make_view<int, 1>(ghost_).data(),
                            ghost_.size());
        size_global_ = gather_scatter_->glb_dof();
    }
}


void PointCloud::adjointHaloExchange(const FieldSet& fieldset, bool on_device) const {
    if (halo_exchange_) {
        for (idx_t f = 0; f < fieldset.size(); ++f) {
            Field& field = const_cast<FieldSet&>(fieldset)[f];
            switch (field.rank()) {
                case 1:
                    dispatch_adjointHaloExchange<1>(field, halo_exchange(), on_device);
                    break;
                case 2:
                    dispatch_adjointHaloExchange<2>(field, halo_exchange(), on_device);
                    break;
                case 3:
                    dispatch_adjointHaloExchange<3>(field, halo_exchange(), on_device);
                    break;
                case 4:
                    dispatch_adjointHaloExchange<4>(field, halo_exchange(), on_device);
                    break;
                default:
                    throw_Exception("Rank not supported", Here());
            }
        }
    }
}

void PointCloud::adjointHaloExchange(const Field& field, bool) const {
    FieldSet fieldset;
    fieldset.add(field);
    adjointHaloExchange(fieldset);
}


}  // namespace detail

PointCloud::PointCloud(const FunctionSpace& functionspace):
    FunctionSpace(functionspace), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const Field& field, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(field,config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const Field& field1, const Field& field2, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(field1, field2, config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const FieldSet& fset, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(fset, config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const std::vector<PointXY>& points, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(points, config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const std::vector<PointXYZ>& points, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(points, config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const std::initializer_list<std::initializer_list<double>>& points, const eckit::Configuration& config):
    FunctionSpace((points.begin()->size() == 2
                       ? new detail::PointCloud{std::vector<PointXY>(points.begin(), points.end()), config}
                       : new detail::PointCloud{std::vector<PointXYZ>(points.begin(), points.end()), config})),
    functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const Grid& grid, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(grid, config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const Grid& grid, const grid::Partitioner& partitioner, const eckit::Configuration& config):
    FunctionSpace(new detail::PointCloud(grid, partitioner,config)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}


}  // namespace functionspace
}  // namespace atlas
