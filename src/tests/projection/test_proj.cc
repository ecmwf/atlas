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

auto check_xy = []( PointXY xy, PointXY ref ) {
    double tolerance_micrometers = 1.e-6;

    auto old = Log::info().precision( 16 );
    Log::info() << "Check " << xy << " = (reference) " << ref << std::endl;
    Log::info().precision( old );

    EXPECT( is_approximately_equal( xy.x(), ref.x(), tolerance_micrometers ) );
    EXPECT( is_approximately_equal( xy.y(), ref.y(), tolerance_micrometers ) );
};

auto check_ll = []( PointLonLat ll, PointLonLat ref ) {
    double tolerance_microdegrees = 1.e-6;

    auto old = Log::info().precision( 16 );
    Log::info() << "Check " << ll << " = (reference) " << ref << std::endl;
    Log::info().precision( old );

    EXPECT( is_approximately_equal( ll.lon(), ref.lon(), tolerance_microdegrees ) );
    EXPECT( is_approximately_equal( ll.lat(), ref.lat(), tolerance_microdegrees ) );
};

//-----------------------------------------------------------------------------

CASE( "test_proj_example_string" ) {
    Projection projection( util::Config( "type", "proj" ).set( "proj", "+proj=utm +zone=32 +datum=WGS84" ) );

    PointLonLat a{12, 55};
    check_ll( a, {12, 55} );

    PointXY b = projection.xy( a );
    check_xy( b, {691875.632137542, 6098907.825129169} );
    check_ll( projection.lonlat( b ), {12, 55} );
}

//-----------------------------------------------------------------------------

CASE( "test_proj_example_EPSG:32632" ) {
    Projection projection( util::Config( "type", "proj" ).set( "proj", "EPSG:32632" ) );

    PointLonLat a{12, 55};
    check_ll( a, {12, 55} );

    PointXY b = projection.xy( a );
    check_xy( b, {691875.632137542, 6098907.825129169} );
    check_ll( projection.lonlat( b ), {12, 55} );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
