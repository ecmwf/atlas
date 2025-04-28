/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

//#include <cstdarg>
//#include <functional>

#include "eckit/utils/MD5.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/detail/Cache.h"

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

class NodeColumnsHaloExchangeCache : public util::Cache<std::string, parallel::HaloExchange>,
                                     public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::HaloExchange>;
    NodeColumnsHaloExchangeCache(): Base("NodeColumnsHaloExchangeCache") {}

public:
    static NodeColumnsHaloExchangeCache& instance() {
        static NodeColumnsHaloExchangeCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const Mesh& mesh, long halo) {
        registerMesh(*mesh.get());
        creator_type creator = std::bind(&NodeColumnsHaloExchangeCache::create, mesh, halo);
        return Base::get_or_create(key(*mesh.get(), halo), creator);
    }
    void onMeshDestruction(mesh::detail::MeshImpl& mesh) override {
        for (long jhalo = 0; jhalo <= mesh::Halo(mesh).size(); ++jhalo) {
            remove(key(mesh, jhalo));
        }
    }

private:
    static Base::key_type key(const mesh::detail::MeshImpl& mesh, long halo) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "],halo[size=" << halo << "]";
        return key.str();
    }

    static value_type* create(const Mesh& mesh, long halo) {
        value_type* value = new value_type();

        std::ostringstream ss;
        ss << "nb_nodes_including_halo[" << halo << "]";
        idx_t nb_nodes(mesh.nodes().size());
        mesh.metadata().get(ss.str(), nb_nodes);

        value->setup(mesh.mpi_comm(), array::make_view<int, 1>(mesh.nodes().partition()).data(),
                     array::make_view<idx_t, 1>(mesh.nodes().remote_index()).data(), REMOTE_IDX_BASE, nb_nodes);

        return value;
    }
};

class NodeColumnsGatherScatterCache : public util::Cache<std::string, parallel::GatherScatter>,
                                      public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::GatherScatter>;
    NodeColumnsGatherScatterCache(): Base("NodeColumnsGatherScatterCache") {}

public:
    static NodeColumnsGatherScatterCache& instance() {
        static NodeColumnsGatherScatterCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const Mesh& mesh) {
        registerMesh(*mesh.get());
        creator_type creator = std::bind(&NodeColumnsGatherScatterCache::create, mesh);
        return Base::get_or_create(key(*mesh.get()), creator);
    }
    void onMeshDestruction(mesh::detail::MeshImpl& mesh) override { remove(key(mesh)); }

private:
    static Base::key_type key(const mesh::detail::MeshImpl& mesh) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "]";
        return key.str();
    }

    static value_type* create(const Mesh& mesh) {
        value_type* value = new value_type();

        mesh::IsGhostNode is_ghost(mesh.nodes());
        std::vector<int> mask(mesh.nodes().size());
        const idx_t npts = mask.size();
        atlas_omp_parallel_for(idx_t n = 0; n < npts; ++n) {
            mask[n] = is_ghost(n) ? 1 : 0;

            // --> This would add periodic west-bc to the gather, but means that
            // global-sums, means, etc are computed wrong
            // if( mask[j] == 1 &&
            // internals::Topology::check(flags(j),internals::Topology::BC) ) {
            //  mask[j] = 0;
            //}
        }

        value->setup(mesh.mpi_comm(),
                     array::make_view<int, 1>(mesh.nodes().partition()).data(),
                     array::make_view<idx_t, 1>(mesh.nodes().remote_index()).data(), REMOTE_IDX_BASE,
                     array::make_view<gidx_t, 1>(mesh.nodes().global_index()).data(), mask.data(), mesh.nodes().size());
        return value;
    }
};

class NodeColumnsChecksumCache : public util::Cache<std::string, parallel::Checksum>,
                                 public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::Checksum>;
    NodeColumnsChecksumCache(): Base("NodeColumnsChecksumCache") {}

public:
    static NodeColumnsChecksumCache& instance() {
        static NodeColumnsChecksumCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const Mesh& mesh) {
        registerMesh(*mesh.get());
        creator_type creator = std::bind(&NodeColumnsChecksumCache::create, mesh);
        return Base::get_or_create(key(*mesh.get()), creator);
    }
    void onMeshDestruction(mesh::detail::MeshImpl& mesh) override { remove(key(mesh)); }

private:
    static Base::key_type key(const mesh::detail::MeshImpl& mesh) {
        std::ostringstream key;
        key << "mesh[address=" << &mesh << "]";
        return key.str();
    }

    static value_type* create(const Mesh& mesh) {
        value_type* value = new value_type();
        util::ObjectHandle<parallel::GatherScatter> gather(
            NodeColumnsGatherScatterCache::instance().get_or_create(mesh));
        value->setup(gather);
        return value;
    }
};

NodeColumns::NodeColumns(Mesh mesh): NodeColumns(mesh, util::NoConfig()) {}

NodeColumns::NodeColumns(Mesh mesh, const eckit::Configuration& config):
    mesh_(mesh), nodes_(mesh_.nodes()), nb_levels_(config.getInt("levels", 0)), nb_nodes_(0) {
    ATLAS_TRACE();
    if (config.has("halo")) {
        halo_ = mesh::Halo(config.getInt("halo"));
    }
    else {
        halo_ = mesh::Halo(mesh);
    }
    mesh::actions::build_nodes_parallel_fields(mesh_);
    mesh::actions::build_periodic_boundaries(mesh_);

    if (halo_.size() > 0) {
        mesh::actions::build_halo(mesh_, halo_.size());
        std::stringstream ss;
        ss << "nb_nodes_including_halo[" << halo_.size() << "]";
        mesh_.metadata().get(ss.str(), nb_nodes_);
    }
    if (!nb_nodes_) {
        std::stringstream ss;
        ss << "nb_nodes_including_halo[" << halo_.size() << "]";
        if (!mesh_.metadata().get(ss.str(), nb_nodes_)) {
            nb_nodes_ = mesh_.nodes().size();
        }
    }

    if (mesh_.grid()) grid_ = mesh_.grid();
}

NodeColumns::~NodeColumns() = default;

std::string NodeColumns::distribution() const {
    return mesh().metadata().getString("distribution");
}

idx_t NodeColumns::nb_nodes() const {
    return nb_nodes_;
}

idx_t NodeColumns::nb_nodes_global() const {
    if (nb_nodes_global_ >= 0) {
        return nb_nodes_global_;
    }
    if (Grid grid = mesh().grid()) {
        nb_nodes_global_ = grid.size();
    }
    else {
        nb_nodes_global_ = gather().glb_dof();
    }
    return nb_nodes_global_;
}

idx_t NodeColumns::config_nb_nodes(const eckit::Configuration& config) const {
    idx_t size = nb_nodes();
    bool global(false);
    if (config.get("global", global)) {
        if (global) {
            idx_t owner(0);
            config.get("owner", owner);
            idx_t _nb_nodes_global = nb_nodes_global();
            idx_t rank = mpi::comm(mpi_comm()).rank();
            size = (rank == owner ? _nb_nodes_global : 0);
        }
    }
    return size;
}

void NodeColumns::set_field_metadata(const eckit::Configuration& config, Field& field) const {
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

    idx_t levels(nb_levels_);
    config.get("levels", levels);
    field.set_levels(levels);

    idx_t variables(0);
    config.get("variables", variables);
    field.set_variables(variables);

    if (config.has("type")) {
      field.metadata().set("type", config.getString("type"));
    }
}

array::DataType NodeColumns::config_datatype(const eckit::Configuration& config) const {
    array::DataType::kind_t kind;
    if (!config.get("datatype", kind)) {
        throw_Exception("datatype missing", Here());
    }
    return array::DataType(kind);
}

std::string NodeColumns::config_name(const eckit::Configuration& config) const {
    std::string name;
    config.get("name", name);
    return name;
}

idx_t NodeColumns::config_levels(const eckit::Configuration& config) const {
    idx_t levels(nb_levels_);
    config.get("levels", levels);
    return levels;
}

array::ArrayShape NodeColumns::config_shape(const eckit::Configuration& config) const {
    array::ArrayShape shape;

    shape.push_back(config_nb_nodes(config));

    idx_t levels(nb_levels_);
    config.get("levels", levels);
    if (levels > 0) {
        shape.push_back(levels);
    }

    idx_t variables(0);
    config.get("variables", variables);
    if (variables > 0) {
        shape.push_back(variables);
    }

    return shape;
}

Field NodeColumns::createField(const eckit::Configuration& config) const {
    Field field = Field(config_name(config), config_datatype(config), config_shape(config));

    set_field_metadata(config, field);

    return field;
}

Field NodeColumns::createField(const Field& other, const eckit::Configuration& config) const {
    return createField(option::name(other.name()) | option::datatype(other.datatype()) | option::levels(other.levels()) |
                       option::variables(other.variables()) | config);
}

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

void NodeColumns::haloExchange(const FieldSet& fieldset, bool on_device) const {
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
    }
}

void NodeColumns::adjointHaloExchange(const FieldSet& fieldset, bool on_device) const {
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

void NodeColumns::haloExchange(const Field& field, bool on_device) const {
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset, on_device);
}

void NodeColumns::adjointHaloExchange(const Field& field, bool) const {
    FieldSet fieldset;
    fieldset.add(field);
    adjointHaloExchange(fieldset);
}

const parallel::HaloExchange& NodeColumns::halo_exchange() const {
    if (halo_exchange_) {
        return *halo_exchange_;
    }
    halo_exchange_ = NodeColumnsHaloExchangeCache::instance().get_or_create(mesh_, halo_.size());
    return *halo_exchange_;
}

void NodeColumns::gather(const FieldSet& local_fieldset, FieldSet& global_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());

    mpi::Scope mpi_scope(mpi_comm());

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

void NodeColumns::gather(const Field& local, Field& global) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields, global_fields);
}
const parallel::GatherScatter& NodeColumns::gather() const {
    if (gather_scatter_) {
        return *gather_scatter_;
    }
    gather_scatter_ = NodeColumnsGatherScatterCache::instance().get_or_create(mesh_);
    return *gather_scatter_;
}
const parallel::GatherScatter& NodeColumns::scatter() const {
    if (gather_scatter_) {
        return *gather_scatter_;
    }
    gather_scatter_ = NodeColumnsGatherScatterCache::instance().get_or_create(mesh_);
    return *gather_scatter_;
}

void NodeColumns::scatter(const FieldSet& global_fieldset, FieldSet& local_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());

    mpi::Scope mpi_scope(mpi_comm());
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

void NodeColumns::scatter(const Field& global, Field& local) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const parallel::Checksum& checksum, const Field& field) {
    auto values = make_leveled_view<const T>(field);
    array::ArrayT<T> surface_field(values.shape(0), values.shape(2));
    auto surface     = array::make_view<T, 2>(surface_field);
    const idx_t npts = values.shape(0);
    atlas_omp_for(idx_t n = 0; n < npts; ++n) {
        for (idx_t j = 0; j < surface.shape(1); ++j) {
            surface(n, j) = 0.;
            for (idx_t l = 0; l < values.shape(1); ++l) {
                surface(n, j) += values(n, l, j);
            }
        }
    }
    return checksum.execute(surface.data(), surface_field.stride(0));
}
}  // namespace

std::string NodeColumns::checksum(const FieldSet& fieldset) const {
    eckit::MD5 md5;
    for (idx_t f = 0; f < fieldset.size(); ++f) {
        const Field& field = fieldset[f];
        if (field.datatype() == array::DataType::kind<int>()) {
            md5 << checksum_3d_field<int>(checksum(), field);
        }
        else if (field.datatype() == array::DataType::kind<long>()) {
            md5 << checksum_3d_field<long>(checksum(), field);
        }
        else if (field.datatype() == array::DataType::kind<float>()) {
            md5 << checksum_3d_field<float>(checksum(), field);
        }
        else if (field.datatype() == array::DataType::kind<double>()) {
            md5 << checksum_3d_field<double>(checksum(), field);
        }
        else {
            throw_Exception("datatype not supported", Here());
        }
    }
    return md5;
}
std::string NodeColumns::checksum(const Field& field) const {
    FieldSet fieldset;
    fieldset.add(field);
    return checksum(fieldset);
}

const parallel::Checksum& NodeColumns::checksum() const {
    if (checksum_) {
        return *checksum_;
    }
    checksum_ = NodeColumnsChecksumCache::instance().get_or_create(mesh_);
    return *checksum_;
}

// std::string NodesFunctionSpace::checksum( const FieldSet& fieldset ) const {
//  const parallel::Checksum& checksum = mesh_.checksum().get(checksum_name());

//  eckit::MD5 md5;
//  for( idx_t f=0; f<fieldset.size(); ++f ) {
//    const Field& field=fieldset[f];
//    if     ( field.datatype() == array::DataType::kind<int>() )
//      md5 << checksum.execute( field.data<int>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<long>() )
//      md5 << checksum.execute( field.data<long>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<float>() )
//      md5 << checksum.execute( field.data<float>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<double>() )
//      md5 << checksum.execute( field.data<double>(), field.stride(0) );
//    else throw_Exception("datatype not supported",Here());
//  }
//  return md5;
//}
// std::string NodesFunctionSpace::checksum( const Field& field ) const {
//  FieldSet fieldset;
//  fieldset.add(field);
//  return checksum(fieldset);
//}

const Grid& NodeColumns::grid() const {
    if (grid_) return grid_;

    const auto& comm = mpi::comm(mpi_comm());
    std::vector<PointXY> points;
    if (comm.size() == 1) {
        const auto xy = atlas::array::make_view<double, 2>(mesh_.nodes().xy());
        for (auto i = 0; i < xy.shape(0); i++) {
            points.push_back({xy(i, 0), xy(i, 1)});
        }
    } else {
        std::vector<int> gidx;
        std::vector<double> x, y;
        const auto gidxView = array::make_view<gidx_t, 1>(global_index());
        const auto ghostView = array::make_view<int, 1>(ghost());
        const auto xy = atlas::array::make_view<double, 2>(mesh_.nodes().xy());
        for (auto i = 0; i < xy.shape(0); i++) {
            if (ghostView(i) == 0) {
                gidx.push_back(gidxView(i));
                x.push_back(xy(i, 0));
                y.push_back(xy(i, 1));
            }
        }
        eckit::mpi::Buffer<int> gidxBuffer(comm.size());
        eckit::mpi::Buffer<double> xBuffer(comm.size());
        eckit::mpi::Buffer<double> yBuffer(comm.size());
        comm.allGatherv(gidx.begin(), gidx.end(), gidxBuffer);
        comm.allGatherv(x.begin(), x.end(), xBuffer);
        comm.allGatherv(y.begin(), y.end(), yBuffer);
        points.reserve(gidxBuffer.buffer.size());
        for (auto i : gidxBuffer.buffer) {
            points[i - 1] = atlas::PointXY{xBuffer.buffer[i - 1], yBuffer.buffer[i - 1]};
        }
    }
    grid_ = UnstructuredGrid(points);
    return grid_;
}

}  // namespace detail

NodeColumns::NodeColumns(): FunctionSpace(), functionspace_(nullptr) {}

NodeColumns::NodeColumns(const FunctionSpace& functionspace):
    FunctionSpace(functionspace), functionspace_(dynamic_cast<const detail::NodeColumns*>(get())) {}

namespace {
detail::NodeColumns* make_functionspace(Mesh mesh, const eckit::Configuration& config) {
    return new detail::NodeColumns(mesh, config);
}
}  // namespace

NodeColumns::NodeColumns(Mesh mesh):
    FunctionSpace(make_functionspace(mesh, util::NoConfig())),
    functionspace_(dynamic_cast<const detail::NodeColumns*>(get())) {}

NodeColumns::NodeColumns(Mesh mesh, const eckit::Configuration& config):
    FunctionSpace(make_functionspace(mesh, config)), functionspace_(dynamic_cast<const detail::NodeColumns*>(get())) {}

idx_t NodeColumns::nb_nodes() const {
    return functionspace_->nb_nodes();
}

idx_t NodeColumns::nb_nodes_global() const {  // All MPI ranks will have same output
    return functionspace_->nb_nodes_global();
}

const Mesh& NodeColumns::mesh() const {
    return functionspace_->mesh();
}

mesh::Nodes& NodeColumns::nodes() const {
    return functionspace_->nodes();
}

// -- Parallelisation aware methods

const mesh::Halo& NodeColumns::halo() const {
    return functionspace_->halo();
}

void NodeColumns::haloExchange(const FieldSet& fieldset, bool on_device) const {
    functionspace_->haloExchange(fieldset, on_device);
}

void NodeColumns::haloExchange(const Field& field, bool on_device) const {
    functionspace_->haloExchange(field, on_device);
}

const parallel::HaloExchange& NodeColumns::halo_exchange() const {
    return functionspace_->halo_exchange();
}

std::string NodeColumns::checksum(const FieldSet& fieldset) const {
    return functionspace_->checksum(fieldset);
}

std::string NodeColumns::checksum(const Field& field) const {
    return functionspace_->checksum(field);
}

const parallel::Checksum& NodeColumns::checksum() const {
    return functionspace_->checksum();
}

}  // namespace functionspace
}  // namespace atlas
