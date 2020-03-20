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

class LonLatPolygon;
class LonLatPolygons;

//------------------------------------------------------------------------------------------------------

/// @brief Implement PolygonCoordinates::contains for a polygon defined in LonLat space.
class LonLatPolygon : public PolygonCoordinates {
private:
    template <typename PointContainer>
    using enable_if_not_polygon = typename std::enable_if<!std::is_base_of<Polygon, PointContainer>::value, int>::type;

public:
    using Vector = LonLatPolygons;

    LonLatPolygon( const Polygon&, const atlas::Field& coordinates, bool removeAlignedPoints = true );
    LonLatPolygon( const PartitionPolygon& );

    template <typename PointContainer, enable_if_not_polygon<PointContainer> = 0>
    LonLatPolygon( const PointContainer& points, bool removeAlignedPoints = true );

    /// @brief Point-in-polygon test based on winding number
    /// @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a Polygon</a>
    /// @param[in] P given point
    /// @return if point is in polygon
    bool contains( const Point2& P ) const override;

private:
    PointLonLat centroid_;
    double inner_radius_squared_{0};
    PointLonLat inner_coordinatesMin_;
    PointLonLat inner_coordinatesMax_;
};


//------------------------------------------------------------------------------------------------------

/// @brief Vector of all polygons with functionality to find partition using a KDTree
class LonLatPolygons : public PolygonCoordinates::Vector {
public:
    LonLatPolygons( const PartitionPolygons& partition_polygons ) {
        reserve( partition_polygons.size() );
        for ( auto& partition_polygon : partition_polygons ) {
            emplace_back( new LonLatPolygon( partition_polygon ) );
        }
    }
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
