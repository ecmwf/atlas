/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Pedro Maciel
/// @author Willem Deconinck
/// @date September 2017

#pragma once

#include <utility>
#include <vector>
#include "atlas/util/Point.h"


namespace atlas {
class Field;
namespace mesh {
namespace detail {
class Polygon;
}
}
}

namespace atlas {
namespace mesh {
namespace detail {

//------------------------------------------------------------------------------------------------------

class PolygonCoordinates {
public:

    // -- Constructors

    PolygonCoordinates(const Polygon&, const atlas::Field& lonlat, bool includesNorthPole, bool includesSouthPole, bool removeAlignedPoints=false);

    PolygonCoordinates(const std::vector<PointLonLat>& points, bool includesNorthPole, bool includesSouthPole);

    // -- Methods

    /*
     * Point-in-partition test based on winding number for a point in a polygon
     * @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a PolygonCoordinates</a>
     * @param[in] points vertex points of a polygon (closed, where poly.front() == poly.back())
     * @param[in] P given point
     * @return if point is in partition
     */
    bool containsPointInLonLatGeometry(const PointLonLat&) const;

    /*
     * Point-in-partition test based on winding number for a point in a polygon
     * @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a PolygonCoordinates</a>
     * @param[in] points vertex points of a polygon (closed, where poly.front() == poly.back())
     * @param[in] P given point
     * @return if point is in partition
     */
    bool containsPointInSphericalGeometry(const PointLonLat&) const;

private:

    // -- Members

    PointLonLat coordinatesMin_;
    PointLonLat coordinatesMax_;
    std::vector< PointLonLat > coordinates_;
    const bool includesNorthPole_;
    const bool includesSouthPole_;

};

//------------------------------------------------------------------------------------------------------

}  // detail
}  // mesh
}  // atlas

