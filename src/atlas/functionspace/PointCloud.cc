/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/functionspace/PointCloud.h"
#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/option/Options.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Metadata.h"

#if ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif

namespace atlas {
namespace functionspace {

namespace detail {

template <>
PointCloud::PointCloud(const std::vector<PointXY>& points) {
    lonlat_     = Field("lonlat", array::make_datatype<double>(), array::make_shape(points.size(), 2));
    auto lonlat = array::make_view<double, 2>(lonlat_);
    for (idx_t j = 0, size = points.size(); j < size; ++j) {
        lonlat(j, 0) = points[j].x();
        lonlat(j, 1) = points[j].y();
    }
}

template <>
PointCloud::PointCloud(const std::vector<PointXYZ>& points) {
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

PointCloud::PointCloud(const Field& lonlat): lonlat_(lonlat) {}

PointCloud::PointCloud(const Field& lonlat, const Field& ghost): lonlat_(lonlat), ghost_(ghost) {}

PointCloud::PointCloud(const FieldSet & flds): lonlat_(flds["lonlat"]),
  ghost_(flds["ghost"]), remote_index_(flds["remote_index"]), partition_(flds["partition"])
{
    ATLAS_ASSERT(ghost_.size() == remote_index_.size());
    ATLAS_ASSERT(ghost_.size() == partition_.size());
    halo_exchange_.reset(new parallel::HaloExchange());
    halo_exchange_->setup(array::make_view<int, 1>( partition_).data(),
                          array::make_view<idx_t, 1>(remote_index_).data(),
                          0,
                          ghost_.size());

    auto lonlat = array::make_view<double, 2>(lonlat_);
    auto ghost = array::make_view<int, 1>(ghost_);
    auto partition = array::make_view<int, 1>(partition_);
    auto remote_index = array::make_view<idx_t, 1>(remote_index_);
    for (idx_t j = 0; j < ghost_.size(); ++j ) {
      std::cout << "pointcloud constructor " << atlas::mpi::rank() << " " << j
                << " lonlat " << lonlat(j, 0) << " " << lonlat(j, 1)  << " "
                << " ghost " << ghost(j) << " "
                << " partition " << partition(j) << " "
                << " remote_index " << remote_index(j) << std::endl;
    }
}

PointCloud::PointCloud(const Grid& grid) {
    lonlat_     = Field("lonlat", array::make_datatype<double>(), array::make_shape(grid.size(), 2));
    auto lonlat = array::make_view<double, 2>(lonlat_);

    idx_t j{0};
    for (auto p : grid.lonlat()) {
        lonlat(j, 0) = p.lon();
        lonlat(j, 1) = p.lat();
        ++j;
    }
}

Field PointCloud::ghost() const {
    if (not ghost_) {
        ghost_ = Field("ghost", array::make_datatype<int>(), array::make_shape(size()));
        array::make_view<int, 1>(ghost_).assign(0);
    }
    return ghost_;
}

array::ArrayShape PointCloud::config_shape(const eckit::Configuration& config) const {
    array::ArrayShape shape;

    shape.emplace_back(size());

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

    auto fldv = array::make_view<int, 2>(field);
    for (idx_t j = 0; j < fldv.shape(0); ++j) {
      std::cout << " PointCloud::haloExchange field before halo " << atlas::mpi::rank() << " "
                << fldv(j, 0)
                << std::endl;
    }

    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset, on_device);

    auto fldv2 = array::make_view<int, 2>(fieldset[field.name()]);
    for (idx_t j = 0; j < fldv2.shape(0); ++j) {
      std::cout << "field after halo " << atlas::mpi::rank() << " "
                << fldv2(j, 0)
                << std::endl;
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

PointCloud::PointCloud(const Field& points):
    FunctionSpace(new detail::PointCloud(points)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const FieldSet& points):
    FunctionSpace(new detail::PointCloud(points)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const std::vector<PointXY>& points):
    FunctionSpace(new detail::PointCloud(points)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const std::vector<PointXYZ>& points):
    FunctionSpace(new detail::PointCloud(points)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const std::initializer_list<std::initializer_list<double>>& points):
    FunctionSpace((points.begin()->size() == 2
                       ? new detail::PointCloud{std::vector<PointXY>(points.begin(), points.end())}
                       : new detail::PointCloud{std::vector<PointXYZ>(points.begin(), points.end())})),
    functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}

PointCloud::PointCloud(const Grid& grid):
    FunctionSpace(new detail::PointCloud(grid)), functionspace_(dynamic_cast<const detail::PointCloud*>(get())) {}


}  // namespace functionspace
}  // namespace atlas
