/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <functional>

#include "eckit/utils/MD5.h"

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/Build2DCellCentres.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/detail/Cache.h"

#include "atlas/field/detail/FieldImpl.h"

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

class CellColumnsHaloExchangeCache : public util::Cache<std::string, parallel::HaloExchange>,
                                     public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::HaloExchange>;
    CellColumnsHaloExchangeCache(): Base("CellColumnsHaloExchangeCache") {}

public:
    static CellColumnsHaloExchangeCache& instance() {
        static CellColumnsHaloExchangeCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const Mesh& mesh) {
        registerMesh(*mesh.get());
        creator_type creator = std::bind(&CellColumnsHaloExchangeCache::create, mesh);
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
        value->setup(mesh.mpi_comm(),
                     array::make_view<int, 1>(mesh.cells().partition()).data(),
                     array::make_view<idx_t, 1>(mesh.cells().remote_index()).data(), REMOTE_IDX_BASE,
                     mesh.cells().size());
        return value;
    }
};

class CellColumnsGatherScatterCache : public util::Cache<std::string, parallel::GatherScatter>,
                                      public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::GatherScatter>;
    CellColumnsGatherScatterCache(): Base("CellColumnsGatherScatterCache") {}

public:
    static CellColumnsGatherScatterCache& instance() {
        static CellColumnsGatherScatterCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const Mesh& mesh) {
        registerMesh(*mesh.get());
        creator_type creator = std::bind(&CellColumnsGatherScatterCache::create, mesh);
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
        value->setup(mesh.mpi_comm(),
                     array::make_view<int, 1>(mesh.cells().partition()).data(),
                     array::make_view<idx_t, 1>(mesh.cells().remote_index()).data(), REMOTE_IDX_BASE,
                     array::make_view<gidx_t, 1>(mesh.cells().global_index()).data(), mesh.cells().size());
        return value;
    }
};

class CellColumnsChecksumCache : public util::Cache<std::string, parallel::Checksum>,
                                 public mesh::detail::MeshObserver {
private:
    using Base = util::Cache<std::string, parallel::Checksum>;
    CellColumnsChecksumCache(): Base("CellColumnsChecksumCache") {}

public:
    static CellColumnsChecksumCache& instance() {
        static CellColumnsChecksumCache inst;
        return inst;
    }
    util::ObjectHandle<value_type> get_or_create(const Mesh& mesh) {
        registerMesh(*mesh.get());
        creator_type creator = std::bind(&CellColumnsChecksumCache::create, mesh);
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
            CellColumnsGatherScatterCache::instance().get_or_create(mesh));
        value->setup(gather);
        return value;
    }
};

void CellColumns::set_field_metadata(const eckit::Configuration& config, Field& field) const {
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
}

idx_t CellColumns::config_size(const eckit::Configuration& config) const {
    idx_t size       = nb_cells();
    bool global(false);
    if (config.get("global", global)) {
        if (global) {
            idx_t owner(0);
            config.get("owner", owner);
            idx_t _nb_cells_global(nb_cells_global());
            const idx_t rank = mpi::comm(mpi_comm()).rank();
            size = (rank == owner ? _nb_cells_global : 0);
        }
    }
    return size;
}

array::DataType CellColumns::config_datatype(const eckit::Configuration& config) const {
    array::DataType::kind_t kind;
    if (!config.get("datatype", kind)) {
        throw_Exception("datatype missing", Here());
    }
    return array::DataType(kind);
}

std::string CellColumns::config_name(const eckit::Configuration& config) const {
    std::string name;
    config.get("name", name);
    return name;
}

idx_t CellColumns::config_levels(const eckit::Configuration& config) const {
    idx_t levels(nb_levels_);
    config.get("levels", levels);
    return levels;
}

array::ArrayShape CellColumns::config_shape(const eckit::Configuration& config) const {
    array::ArrayShape shape;

    shape.push_back(config_size(config));

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

CellColumns::CellColumns(const Mesh& mesh, const eckit::Configuration& config):
    mesh_(mesh), cells_(mesh_.cells()), nb_levels_(config.getInt("levels", 0)), nb_cells_(0) {
    ATLAS_TRACE();
    if (config.has("halo")) {
        halo_ = mesh::Halo(config.getInt("halo"));
    }
    else {
        halo_ = mesh::Halo(mesh_);
    }

    auto get_nb_cells_from_metadata = [&]() {
        idx_t nb_cells{0};
        for (idx_t t = 0; t < cells_.nb_types(); ++t) {
            std::stringstream ss;
            ss << "nb_cells_including_halo[" << t << "][" << halo_.size() << "]";
            idx_t nb_cells_for_this_type{0};
            mesh_.metadata().get(ss.str(), nb_cells_for_this_type);
            nb_cells += nb_cells_for_this_type;
        }
        return nb_cells;
    };

    mesh::actions::build_nodes_parallel_fields(mesh_);
    mesh::actions::build_cells_parallel_fields(mesh_);
    mesh::actions::build_periodic_boundaries(mesh_);


    mesh::actions::build_halo(mesh_, halo_.size());
    nb_cells_ = get_nb_cells_from_metadata();

    if (!nb_cells_) {
        nb_cells_ = mesh.cells().size();
    }
    ATLAS_ASSERT(nb_cells_);

    if (mesh_.grid()) {
        grid_ = mesh_.grid();
    }
}

CellColumns::~CellColumns() = default;

size_t CellColumns::footprint() const {
    size_t size = sizeof(*this);
    // TODO
    return size;
}

std::string CellColumns::distribution() const {
    return mesh().metadata().getString("distribution");
}

idx_t CellColumns::nb_cells() const {
    return nb_cells_;
}

idx_t CellColumns::nb_cells_global() const {
    if (nb_cells_global_ >= 0) {
        return nb_cells_global_;
    }
    nb_cells_global_ = gather().glb_dof();
    return nb_cells_global_;
}

idx_t CellColumns::levels() const {
    return nb_levels_;
}

Field CellColumns::createField(const eckit::Configuration& options) const {
    Field field(config_name(options), config_datatype(options), config_shape(options));
    set_field_metadata(options, field);
    return field;
}

Field CellColumns::createField(const Field& other, const eckit::Configuration& config) const {
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
}  // namespace

void CellColumns::haloExchange(const FieldSet& fieldset, bool on_device) const {
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
void CellColumns::haloExchange(const Field& field, bool on_device) const {
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset, on_device);
}
const parallel::HaloExchange& CellColumns::halo_exchange() const {
    if (halo_exchange_) {
        return *halo_exchange_;
    }
    halo_exchange_ = CellColumnsHaloExchangeCache::instance().get_or_create(mesh_);
    return *halo_exchange_;
}

void CellColumns::gather(const FieldSet& local_fieldset, FieldSet& global_fieldset) const {
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

void CellColumns::gather(const Field& local, Field& global) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields, global_fields);
}
const parallel::GatherScatter& CellColumns::gather() const {
    if (gather_scatter_) {
        return *gather_scatter_;
    }
    gather_scatter_ = CellColumnsGatherScatterCache::instance().get_or_create(mesh_);
    return *gather_scatter_;
}
const parallel::GatherScatter& CellColumns::scatter() const {
    if (gather_scatter_) {
        return *gather_scatter_;
    }
    gather_scatter_ = CellColumnsGatherScatterCache::instance().get_or_create(mesh_);
    return *gather_scatter_;
}

void CellColumns::scatter(const FieldSet& global_fieldset, FieldSet& local_fieldset) const {
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
void CellColumns::scatter(const Field& global, Field& local) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const parallel::Checksum& checksum, const Field& field) {
    auto values = array::make_view<const T, 3>(field);
    array::ArrayT<T> surface_field(field.shape(0), field.shape(2));
    auto surface = array::make_view<T, 2>(surface_field);
    for (idx_t n = 0; n < values.shape(0); ++n) {
        for (idx_t j = 0; j < surface.shape(1); ++j) {
            surface(n, j) = 0.;
            for (idx_t l = 0; l < values.shape(1); ++l) {
                surface(n, j) += values(n, l, j);
            }
        }
    }
    return checksum.execute(surface.data(), surface_field.stride(0));
}
template <typename T>
std::string checksum_2d_field(const parallel::Checksum& checksum, const Field& field) {
    auto values = array::make_view<T, 2>(field);
    return checksum.execute(values.data(), field.stride(0));
}
template <typename T>
std::string checksum_1d_field(const parallel::Checksum& checksum, const Field& field) {
    auto values = array::make_view<T, 1>(field);
    return checksum.execute(values.data(), 1);
}

}  // namespace

std::string CellColumns::checksum(const FieldSet& fieldset) const {
    eckit::MD5 md5;
    for (idx_t f = 0; f < fieldset.size(); ++f) {
        const Field& field = fieldset[f];
        std::string field_checksum;
        if (field.datatype() == array::DataType::kind<int>()) {
            if (field.rank()==3) {
                field_checksum = checksum_3d_field<int>(checksum(), field);
            }
            else if (field.rank()==2) {
                field_checksum = checksum_2d_field<int>(checksum(), field);
            }
            else if (field.rank()==1) {
                field_checksum = checksum_1d_field<int>(checksum(), field);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        else if (field.datatype() == array::DataType::kind<long>()) {
            if (field.rank()==3) {
                field_checksum = checksum_3d_field<long>(checksum(), field);
            }
            else if (field.rank()==2) {
                field_checksum = checksum_2d_field<long>(checksum(), field);
            }
            else if (field.rank()==1) {
                field_checksum = checksum_1d_field<long>(checksum(), field);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        else if (field.datatype() == array::DataType::kind<float>()) {
            if (field.rank()==3) {
                field_checksum = checksum_3d_field<float>(checksum(), field);
            }
            else if (field.rank()==2) {
                field_checksum = checksum_2d_field<float>(checksum(), field);
            }
            else if (field.rank()==1) {
                field_checksum = checksum_1d_field<float>(checksum(), field);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        else if (field.datatype() == array::DataType::kind<double>()) {
            if (field.rank()==3) {
                field_checksum = checksum_3d_field<double>(checksum(), field);
            }
            else if (field.rank()==2) {
                field_checksum = checksum_2d_field<double>(checksum(), field);
            }
            else if (field.rank()==1) {
                field_checksum = checksum_1d_field<double>(checksum(), field);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        else {
            throw_Exception("datatype not supported", Here());
        }
        if (fieldset.size() == 1) {
            return field_checksum;
        }
        else {
            md5 << field_checksum;
        }
    }
    return md5;
}
std::string CellColumns::checksum(const Field& field) const {
    FieldSet fieldset;
    fieldset.add(field);
    return checksum(fieldset);
}

const parallel::Checksum& CellColumns::checksum() const {
    if (checksum_) {
        return *checksum_;
    }
    checksum_ = CellColumnsChecksumCache::instance().get_or_create(mesh_);
    return *checksum_;
}

const Grid& CellColumns::grid() const {
    if (grid_) {
        return grid_;
    }

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

Field CellColumns::lonlat() const {
    if (!mesh_.cells().has_field("lonlat")) {
        mesh::actions::Build2DCellCentres("lonlat")(const_cast<Mesh&>(mesh_));
    }
    return mesh_.cells().field("lonlat");
}

Field CellColumns::remote_index() const {
    return mesh_.cells().remote_index();
}

Field CellColumns::global_index() const {
    return mesh_.cells().global_index();
}

Field CellColumns::ghost() const {
    if (mesh_.cells().has_field("ghost")) {
       return mesh_.cells().field("ghost");
    }
    return mesh_.cells().halo();
}

Field CellColumns::partition() const {
    return mesh_.cells().partition();
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
extern "C" {
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

CellColumns* atlas__fs__CellColumns__new(Mesh::Implementation* mesh, const eckit::Configuration* config) {
    ATLAS_ASSERT(mesh != nullptr);
    Mesh m(mesh);
    return new CellColumns(m, *config);
}

//------------------------------------------------------------------------------

void atlas__fs__CellColumns__delete(CellColumns* This) {
    ATLAS_ASSERT(This != nullptr);
    delete (This);
}

//------------------------------------------------------------------------------

int atlas__fs__CellColumns__nb_cells(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr);
    return This->nb_cells();
}

//------------------------------------------------------------------------------

Mesh::Implementation* atlas__fs__CellColumns__mesh(CellColumns* This) {
    ATLAS_ASSERT(This != nullptr);
    return This->mesh().get();
}

//------------------------------------------------------------------------------

mesh::Cells* atlas__fs__CellColumns__cells(CellColumns* This) {
    ATLAS_ASSERT(This != nullptr);
    return &This->cells();
}

//------------------------------------------------------------------------------

using field::FieldImpl;
using field::FieldSetImpl;

field::FieldImpl* atlas__fs__CellColumns__create_field(const CellColumns* This, const eckit::Configuration* options) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(options);
    FieldImpl* field;
    {
        Field f = This->createField(*options);
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

//------------------------------------------------------------------------------

field::FieldImpl* atlas__fs__CellColumns__create_field_template(const CellColumns* This,
                                                                const field::FieldImpl* field_template,
                                                                const eckit::Configuration* options) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(options);
    FieldImpl* field;
    {
        Field f = This->createField(Field(field_template), *options);
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__halo_exchange_fieldset(const CellColumns* This, field::FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(fieldset != nullptr);
    FieldSet f(fieldset);
    This->haloExchange(f);
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__halo_exchange_field(const CellColumns* This, field::FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(field != nullptr);
    Field f(field);
    This->haloExchange(f);
}

// -----------------------------------------------------------------------------------

const parallel::HaloExchange* atlas__fs__CellColumns__get_halo_exchange(const CellColumns* This) {
    ATLAS_ASSERT(This != nullptr);
    return &This->halo_exchange();
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__gather_fieldset(const CellColumns* This, const field::FieldSetImpl* local,
                                             field::FieldSetImpl* global) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(local);
    ATLAS_ASSERT(global);
    const FieldSet l(local);
    FieldSet g(global);
    This->gather(l, g);
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__gather_field(const CellColumns* This, const field::FieldImpl* local,
                                          field::FieldImpl* global) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(local);
    ATLAS_ASSERT(global);
    const Field l(local);
    Field g(global);
    This->gather(l, g);
}

// -----------------------------------------------------------------------------------

const parallel::GatherScatter* atlas__fs__CellColumns__get_gather(const CellColumns* This) {
    ATLAS_ASSERT(This);
    return &This->gather();
}

// -----------------------------------------------------------------------------------

const parallel::GatherScatter* atlas__fs__CellColumns__get_scatter(const CellColumns* This) {
    ATLAS_ASSERT(This);
    return &This->scatter();
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__scatter_fieldset(const CellColumns* This, const field::FieldSetImpl* global,
                                              field::FieldSetImpl* local) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(local);
    ATLAS_ASSERT(global);
    const FieldSet g(global);
    FieldSet l(local);
    This->scatter(g, l);
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__scatter_field(const CellColumns* This, const field::FieldImpl* global,
                                           field::FieldImpl* local) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(global);
    ATLAS_ASSERT(local);
    const Field g(global);
    Field l(local);
    This->scatter(g, l);
}

// -----------------------------------------------------------------------------------

const parallel::Checksum* atlas__fs__CellColumns__get_checksum(const CellColumns* This) {
    ATLAS_ASSERT(This);
    return &This->checksum();
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__checksum_fieldset(const CellColumns* This, const field::FieldSetImpl* fieldset,
                                               char*& checksum, int& size, int& allocated) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(fieldset);
    std::string checksum_str(This->checksum(fieldset));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}

// -----------------------------------------------------------------------------------

void atlas__fs__CellColumns__checksum_field(const CellColumns* This, const field::FieldImpl* field, char*& checksum,
                                            int& size, int& allocated) {
    ATLAS_ASSERT(This);
    ATLAS_ASSERT(field);
    std::string checksum_str(This->checksum(field));
    size      = static_cast<int>(checksum_str.size());
    checksum  = new char[size + 1];
    allocated = true;
    std::strncpy(checksum, checksum_str.c_str(), size + 1);
}
}

// -----------------------------------------------------------------------------------

}  // namespace detail

// -----------------------------------------------------------------------------------

CellColumns::CellColumns(): FunctionSpace(), functionspace_(nullptr) {}

CellColumns::CellColumns(const FunctionSpace& functionspace):
    FunctionSpace(functionspace), functionspace_(dynamic_cast<const detail::CellColumns*>(get())) {}

CellColumns::CellColumns(const Mesh& mesh, const eckit::Configuration& config):
    FunctionSpace(new detail::CellColumns(mesh, config)),
    functionspace_(dynamic_cast<const detail::CellColumns*>(get())) {}

CellColumns::CellColumns(const Mesh& mesh):
    FunctionSpace(new detail::CellColumns(mesh)), functionspace_(dynamic_cast<const detail::CellColumns*>(get())) {}

idx_t CellColumns::nb_cells() const {
    return functionspace_->nb_cells();
}

idx_t CellColumns::nb_cells_global() const {  // Only on MPI rank 0, will this be different from 0
    return functionspace_->nb_cells_global();
}

idx_t CellColumns::levels() const {
    return functionspace_->levels();
}

const Mesh& CellColumns::mesh() const {
    return functionspace_->mesh();
}

const mesh::HybridElements& CellColumns::cells() const {
    return functionspace_->cells();
}

const parallel::HaloExchange& CellColumns::halo_exchange() const {
    return functionspace_->halo_exchange();
}

std::string CellColumns::checksum(const FieldSet& fieldset) const {
    return functionspace_->checksum(fieldset);
}

std::string CellColumns::checksum(const Field& field) const {
    return functionspace_->checksum(field);
}

const parallel::Checksum& CellColumns::checksum() const {
    return functionspace_->checksum();
}

const mesh::Halo& CellColumns::halo() const {
    return functionspace_->halo();
}

}  // namespace functionspace
}  // namespace atlas
