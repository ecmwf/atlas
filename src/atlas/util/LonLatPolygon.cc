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
#include <iostream>
#include <limits>

#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatPolygon.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

namespace {

double cross_product_analog( const PointLonLat& A, const PointLonLat& B, const PointLonLat& C ) {
    return ( A.lon() - C.lon() ) * ( B.lat() - C.lat() ) - ( A.lat() - C.lat() ) * ( B.lon() - C.lon() );
}

}  // namespace

//------------------------------------------------------------------------------------------------------

LonLatPolygon::LonLatPolygon( const Polygon& poly, const atlas::Field& lonlat, bool removeAlignedPoints ) :
    PolygonCoordinates( poly, lonlat, removeAlignedPoints ) {}

LonLatPolygon::LonLatPolygon( const std::vector<PointLonLat>& points ) : PolygonCoordinates( points ) {}

bool LonLatPolygon::contains( const PointLonLat& P ) const {
    ATLAS_ASSERT( coordinates_.size() >= 2 );

    // check first bounding box
    if ( coordinatesMax_.lon() <= P.lon() || P.lon() < coordinatesMin_.lon() || coordinatesMax_.lat() <= P.lat() ||
         P.lat() < coordinatesMin_.lat() ) {
        return false;
    }

    // winding number
    int wn = 0;

    // loop on polygon edges
    for ( size_t i = 1; i < coordinates_.size(); ++i ) {
        const PointLonLat& A = coordinates_[i - 1];
        const PointLonLat& B = coordinates_[i];

        // check point-edge side and direction, using 2D-analog cross-product;
        // tests if P is left|on|right of a directed A-B infinite line, by
        // intersecting either:
        // - "up" on upward crossing & P left of edge, or
        // - "down" on downward crossing & P right of edge
        const bool APB = ( A.lat() <= P.lat() && P.lat() < B.lat() );
        const bool BPA = ( B.lat() <= P.lat() && P.lat() < A.lat() );

        if ( APB != BPA ) {
            const double side = cross_product_analog( P, A, B );
            if ( APB && side > 0 ) {
                ++wn;
            }
            else if ( BPA && side < 0 ) {
                --wn;
            }
        }
    }

    // wn == 0 only when P is outside
    return wn != 0;
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
