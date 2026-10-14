/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/field/builders/GridBoxArea.h"

#include <algorithm>
#include <cmath>

#include "atlas/array.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/util/Geometry.h"

namespace atlas {
namespace field {
namespace builders {

namespace {
class StructuredColumnsAreaComputer {
public:
    StructuredColumnsAreaComputer(const functionspace::StructuredColumns& fs,
                                  double radius = geometry::Earth().radius());

    double operator()(idx_t point_index) const;
    void operator()(double area[], std::size_t size) const;

private:
    double sin_latitude_strip_width(double lat_a, double lat_b) const;

private:
    const functionspace::StructuredColumns& fs_;
    double radius_;
    array::ArrayView<const idx_t, 1> index_i_;
    array::ArrayView<const idx_t, 1> index_j_;
    bool projection_in_degrees_;
};

StructuredColumnsAreaComputer::StructuredColumnsAreaComputer(const functionspace::StructuredColumns& fs,
                                                             double radius):
    fs_(fs),
    radius_(radius),
    index_i_(array::make_view<const idx_t, 1>(fs.index_i())),
    index_j_(array::make_view<const idx_t, 1>(fs.index_j())) {
    auto proj_type = fs_.projection().type();
    if (proj_type == "rotated_lonlat" || proj_type == "lonlat") {
        projection_in_degrees_ = true;
    }
    else if (fs_.projection().units() == "meters") {
        projection_in_degrees_ = false;
    }
    else {
        ATLAS_DEBUG("Projection type = " << proj_type);
        throw_Exception("StructuredColumnsAreaComputer only works for (rotated) lonlat grids or projections in meters",
                        Here());
    }
}

double StructuredColumnsAreaComputer::operator()(idx_t point_index) const {
    idx_t i          = index_i_(point_index);
    idx_t j          = index_j_(point_index);
    PointXY xy       = fs_.compute_xy(i, j);
    double y         = xy.y();
    double yN        = 0.5 * (y + fs_.compute_xy(i, j - 1).y());
    double yS        = 0.5 * (y + fs_.compute_xy(i, j + 1).y());
    idx_t j_internal = fs_.compute_j(j);
    const double dx  = fs_.grid().dx(j_internal);
    if (projection_in_degrees_) {
        constexpr double deg2rad = M_PI / 180.0;
        return radius_ * radius_ * dx * deg2rad * sin_latitude_strip_width(yN, yS);
    }
    return dx * std::abs(yN - yS);
}

void StructuredColumnsAreaComputer::operator()(double area[], std::size_t size) const {
    ATLAS_ASSERT(size == fs_.size());
    for (idx_t j = fs_.j_begin_halo(); j < fs_.j_end_halo(); ++j) {
        const idx_t i_begin = fs_.i_begin_halo(j);
        const idx_t i_end   = fs_.i_end_halo(j);
        if (i_begin == i_end) {
            continue;
        }

        const idx_t first_row_point       = fs_.index(i_begin, j);
        const double first_row_point_area = operator()(first_row_point);

        for (idx_t i = i_begin; i < i_end; ++i) {
            area[fs_.index(i, j)] = first_row_point_area;
        }
    }
}

double StructuredColumnsAreaComputer::sin_latitude_strip_width(double lat_a, double lat_b) const {
    constexpr double deg2rad = M_PI / 180.0;
    const double lat_max     = std::max(lat_a, lat_b);
    const double lat_min     = std::min(lat_a, lat_b);
    if (lat_max > 90. && lat_min < 90.) {
        return std::abs(std::sin(lat_max * deg2rad) - std::sin(90. * deg2rad)) +
               std::abs(std::sin(lat_min * deg2rad) - std::sin(90. * deg2rad));
    }
    if (lat_min < -90. && lat_max > -90.) {
        return std::abs(std::sin(lat_max * deg2rad) - std::sin(-90. * deg2rad)) +
               std::abs(std::sin(lat_min * deg2rad) - std::sin(-90. * deg2rad));
    }
    return std::abs(std::sin(lat_max * deg2rad) - std::sin(lat_min * deg2rad));
}

}  // namespace

static FieldBuilderFactoryBuilder<GridBoxArea> register_grid_box_area("grid-box-area");

GridBoxArea::GridBoxArea(const util::Config& config, const util::ObjectHandleBase& object) {
    const auto* fs = dynamic_cast<const FunctionSpace*>(&object);
    if (not fs) {
        throw_Exception("GridBoxArea requires a FunctionSpace object handle", Here());
    }
    fs_ = functionspace::StructuredColumns(*fs);
    if (not fs_) {
        throw_Exception("GridBoxArea only works for StructuredColumns function spaces", Here());
    }
    field_name_ = config.getString("field_name", "area");
    radius_     = geometry::Earth().radius();
    config.get("radius", radius_);
    if (config.has("geometry")) {
        const std::string geometry_name = config.getString("geometry");
        radius_                         = Geometry(geometry_name).radius();
    }
}

Field GridBoxArea::operator()() const {
    Field area = fs_.createField<double>(option::name(field_name_) | option::levels(0));
    StructuredColumnsAreaComputer area_computer(fs_, radius_);
    area_computer(array::make_view<double, 1>(area).data(), area.size());
    return area;
}

}  // namespace builders
}  // namespace field
}  // namespace atlas