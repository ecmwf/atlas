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

#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/SphericalPolygon.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

SphericalPolygon::SphericalPolygon( const Polygon& poly, const atlas::Field& lonlat ) :
    PolygonCoordinates( poly, lonlat, false ) {}

SphericalPolygon::SphericalPolygon( const std::vector<PointLonLat>& points ) : PolygonCoordinates( points ) {}

bool SphericalPolygon::contains( const PointLonLat& P ) const {
    using eckit::types::is_approximately_equal;

    ATLAS_ASSERT( coordinates_.size() >= 2 );

    // check first bounding box
    if ( coordinatesMax_.lon() <= P.lon() || P.lon() < coordinatesMin_.lon() ) {
        return false;
    }

    // winding number
    int wn = 0;

    // loop on polygon edges
    for ( size_t i = 1; i < coordinates_.size(); ++i ) {
        const PointLonLat& A = coordinates_[i - 1];
        const PointLonLat& B = coordinates_[i];

        // test if P is on/above/below of a great circle containing A,B
        const bool APB = ( A.lon() <= P.lon() && P.lon() < B.lon() );
        const bool BPA = ( B.lon() <= P.lon() && P.lon() < A.lon() );

        if ( APB != BPA ) {
            const double lat = [&] {
                if ( is_approximately_equal( A.lat(), B.lat() ) &&
                     is_approximately_equal( std::abs( A.lat() ), 90. ) ) {
                    return A.lat();
                }
                else {
                    return util::Earth::greatCircleLatitudeGivenLongitude( A, B, P.lon() );
                }
            }();

            ATLAS_ASSERT( !std::isnan( lat ) );
            if ( is_approximately_equal( P.lat(), lat ) ) {
                return true;
            }

            wn += ( P.lat() > lat ? -1 : 1 ) * ( APB ? -1 : 1 );
        }
    }

    // wn == 0 only when P is outside
    return wn != 0;
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
