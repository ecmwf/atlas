/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <type_traits>

#include "atlas/util/Polygon.h"

namespace atlas {
namespace util {

class PartitionPolygon;
class PolygonXY;
class ListPolygonXY;

//------------------------------------------------------------------------------------------------------

/// @brief Implement PolygonCoordinates::contains for a polygon defined in XY space.
class PolygonXY : public PolygonCoordinates {
public:
    using Vector = ListPolygonXY;

    PolygonXY(const PartitionPolygon&);

    /// @brief Point-in-polygon test based on winding number
    /// @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a Polygon</a>
    /// @param[in] P given point
    /// @return if point (x,y) is in polygon
    bool contains(const Point2& Pxy) const override;

private:
    Point2 centroid_;
    double inner_radius_squared_{0};
    Point2 inner_coordinatesMin_;
    Point2 inner_coordinatesMax_;
};


//------------------------------------------------------------------------------------------------------

/// @brief Vector of all polygons with functionality to find partition using a KDTree
class ListPolygonXY : public PolygonCoordinates::Vector {
public:
    ListPolygonXY(const PartitionPolygons& partition_polygons) {
        reserve(partition_polygons.size());
        for (auto& partition_polygon : partition_polygons) {
            emplace_back(new PolygonXY(partition_polygon));
        }
    }
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
