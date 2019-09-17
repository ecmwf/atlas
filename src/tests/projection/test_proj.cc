/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/projection.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

template <typename POINT>
void check( const POINT& a, const POINT& b ) {
    static constexpr double eps = 1.e-6;

    auto old = Log::info().precision( 16 );
    Log::info() << "Check " << a << " = " << b << std::endl;
    Log::info().precision( old );

    EXPECT( is_approximately_equal( a[0], b[0], eps ) );
    EXPECT( is_approximately_equal( a[1], b[1], eps ) );
}

//-----------------------------------------------------------------------------

CASE( "test_proj" ) {
    PointLonLat a{12, 55};
    check( a, {12, 55} );

    for ( auto proj : {"+proj=utm +zone=32 +datum=WGS84", "EPSG:32632"} ) {
        Projection projection( util::Config( "type", "proj" ).set( "proj", proj ) );

        PointXY b = projection.xy( a );
        check( b, {691875.632137542, 6098907.825129169} );
        check( projection.lonlat( b ), a );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
