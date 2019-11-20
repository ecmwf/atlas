/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <limits>
#include <vector>

#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Earth;

namespace atlas {
namespace test {


CASE( "test PointLonLat normalisation" ) {
    // latitude at Valparaíso-Shanghai mid-point
    PointLonLat p;

    p = PointLonLat( -71.6, -33. );
    p.normalise();
    EXPECT( is_approximately_equal( p.lon(), -71.6 + 360. ) );
    p.normalise( -180., 180. );
    EXPECT( is_approximately_equal( p.lon(), -71.6 ) );

    p = PointLonLat( 121.8, 31.4 );
    p.normalise( -180., 180. );
    EXPECT( is_approximately_equal( p.lon(), 121.8 ) );

    p = PointLonLat( 181., 31.4 );
    p.normalise( -180., 180. );
    EXPECT( is_approximately_equal( p.lon(), 181. - 360. ) );

    p = PointLonLat( 180., 31.4 );
    p.normalise( -180., 180. );
    EXPECT( is_approximately_equal( p.lon(), 180. ) );

    p = PointLonLat( 180., 31.4 );
    p.normalise( -180 );
    EXPECT( is_approximately_equal( p.lon(), -180. ) );

    p = PointLonLat( -180., 31.4 );
    p.normalise( -180., 180. );
    EXPECT( is_approximately_equal( p.lon(), -180. ) );

    p = PointLonLat( -180., 31.4 );
    p.normalise( -180. );
    EXPECT( is_approximately_equal( p.lon(), -180. ) );
}

CASE( "test output vector<PointXY>" ) {
    // clang-format off
    std::vector<PointXY> points = {
        {0.,1.},
        {0.,2.},
        {0.,3.},
        {1.,1.},
        {1.,2.},
        {1.,3.},
        {2.,1.},
        {3.,2.},
        {4.,3.},
    };
    // clang-format on
    Log::info() << "points = " << points << std::endl;
}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
