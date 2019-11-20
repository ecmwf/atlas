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


//------------------------------------------------------------------------------------------------------

class LonLatPolygon : public PolygonCoordinates {
private:
    template <typename PointContainer>
    using enable_if_not_polygon = typename std::enable_if<!std::is_base_of<Polygon, PointContainer>::value, int>::type;

public:
    // -- Constructors

    LonLatPolygon( const Polygon&, const atlas::Field& coordinates, bool removeAlignedPoints = true );
    LonLatPolygon( const PartitionPolygon& );

    template <typename PointContainer, enable_if_not_polygon<PointContainer> = 0>
    LonLatPolygon( const PointContainer& points, bool removeAlignedPoints = true );

    // -- Overridden methods

    /*
   * Point-in-polygon test based on winding number
   * @note reference <a
   * href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point
   * in a Polygon</a>
   * @param[in] P given point
   * @return if point is in polygon
   */
    bool contains( const Point2& P ) const;

private:
    PointLonLat centroid_;
    double inner_radius_squared_{0};
    PointLonLat inner_coordinatesMin_;
    PointLonLat inner_coordinatesMax_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
