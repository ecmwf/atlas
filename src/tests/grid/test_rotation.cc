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

#include "eckit/types/FloatCompare.h"

#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Rotation.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Config;
using atlas::util::Rotation;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

constexpr double eps() {
    return 1.e-5;
}
const double d2r = atlas::util::Constants::degreesToRadians();
const double r2d = atlas::util::Constants::radiansToDegrees();

#define EXPECT_EQUIVALENT( p1, p2 )                                                                \
    if ( std::abs( p2.lat() ) == 90. )                                                             \
        EXPECT( eckit::types::is_approximately_equal( p1.lat(), p2.lat(), eps() ) );               \
    else {                                                                                         \
        EXPECT( eckit::types::is_approximately_equal( 10. + std::cos( p1.lon() * d2r ),            \
                                                      10. + std::cos( p2.lon() * d2r ), eps() ) ); \
        EXPECT( eckit::types::is_approximately_equal( 10. + std::sin( p1.lat() * d2r ),            \
                                                      10. + std::sin( p2.lat() * d2r ), eps() ) ); \
    }

//-----------------------------------------------------------------------------

class MagicsRotation {
    // For reference, this what magics uses, it appears as if it originated from
    // fortran code
    // Strangely the definition of rotate and unrotate are switched.

public:
    MagicsRotation( const PointLonLat& south_pole ) : south_pole_( south_pole ) {}

    PointLonLat rotate( const PointLonLat& point ) const {
        return magics_unrotate( point );  /// Switch meaning !!!
    }

    PointLonLat unrotate( const PointLonLat& point ) const {
        return magics_rotate( point );  /// Swich meaning !!!
    }

private:
    PointLonLat south_pole_;

    PointLonLat magics_rotate( const PointLonLat& point ) const {
        double lat_y = point.lat();
        double lon_x = point.lon();

        double sin_south_pole_lat = std::sin( d2r * ( south_pole_.lat() + 90. ) );
        double cos_south_pole_lat = std::cos( d2r * ( south_pole_.lat() + 90. ) );

        double ZXMXC           = d2r * ( lon_x - south_pole_.lon() );
        double sin_lon_decr_sp = std::sin( ZXMXC );
        double cos_lon_decr_sp = std::cos( ZXMXC );
        double sin_lat         = std::sin( d2r * lat_y );
        double cos_lat         = std::cos( d2r * lat_y );
        double ZSYROT          = cos_south_pole_lat * sin_lat - sin_south_pole_lat * cos_lat * cos_lon_decr_sp;
        ZSYROT                 = std::max( std::min( ZSYROT, +1.0 ), -1.0 );

        double PYROT = std::asin( ZSYROT ) * r2d;

        double ZCYROT = std::cos( PYROT * d2r );
        double ZCXROT = ( cos_south_pole_lat * cos_lat * cos_lon_decr_sp + sin_south_pole_lat * sin_lat ) / ZCYROT;
        ZCXROT        = std::max( std::min( ZCXROT, +1.0 ), -1.0 );
        double ZSXROT = cos_lat * sin_lon_decr_sp / ZCYROT;

        double PXROT = std::acos( ZCXROT ) * r2d;

        if ( ZSXROT < 0.0 ) PXROT = -PXROT;

        return PointLonLat( PXROT, PYROT );
    }

    PointLonLat magics_unrotate( const PointLonLat& point ) const {
        double lat_y = point.lat();
        double lon_x = point.lon();

        double sin_south_pole_lat = std::sin( d2r * ( south_pole_.lat() + 90. ) );
        double cos_south_pole_lat = std::cos( d2r * ( south_pole_.lat() + 90. ) );
        double cos_lon            = std::cos( d2r * lon_x );
        double sin_lat            = std::sin( d2r * lat_y );
        double cos_lat            = std::cos( d2r * lat_y );
        double ZSYREG             = cos_south_pole_lat * sin_lat + sin_south_pole_lat * cos_lat * cos_lon;
        ZSYREG                    = std::max( std::min( ZSYREG, +1.0 ), -1.0 );
        double PYREG              = std::asin( ZSYREG ) * r2d;
        double ZCYREG             = std::cos( PYREG * d2r );
        double ZCXMXC             = ( cos_south_pole_lat * cos_lat * cos_lon - sin_south_pole_lat * sin_lat ) / ZCYREG;
        ZCXMXC                    = std::max( std::min( ZCXMXC, +1.0 ), -1.0 );
        double ZSXMXC             = cos_lat * sin_lat / ZCYREG;
        double ZXMXC              = std::acos( ZCXMXC ) * r2d;
        if ( ZSXMXC < 0.0 ) ZXMXC = -ZXMXC;
        double PXREG = ZXMXC + south_pole_.lon();

        return PointLonLat( PXREG, PYREG );
    }
};

//-----------------------------------------------------------------------------

CASE( "test_rotation" ) {
    Config config;
    config.set( "north_pole", std::vector<double>{-176, 40} );
    Rotation rotation( config );
    MagicsRotation magics( rotation.southPole() );
    Log::info() << rotation << std::endl;

    EXPECT( rotation.rotated() );

    PointLonLat p, r;

    p = {0., 90.};
    r = {-176., 40.};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( magics.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
    EXPECT_EQUIVALENT( magics.unrotate( r ), p );

    p = {0., 0.};
    r = {-176., -50.};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( magics.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
    EXPECT_EQUIVALENT( magics.unrotate( r ), p );

    p = {-180., 45.};
    r = {-176., 85.};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( magics.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
    EXPECT_EQUIVALENT( magics.unrotate( r ), p );
}

CASE( "test_no_rotation" ) {
    Config config;
    Rotation rotation( config );
    MagicsRotation magics( rotation.southPole() );

    Log::info() << rotation << std::endl;

    EXPECT( not rotation.rotated() );

    PointLonLat p, r;

    p = {0., 90.};
    r = p;
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( magics.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
    EXPECT_EQUIVALENT( magics.unrotate( r ), p );

    p = {0., 0.};
    r = p;
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( magics.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
    EXPECT_EQUIVALENT( magics.unrotate( r ), p );

    p = {-180., 45.};
    r = p;
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( magics.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
    EXPECT_EQUIVALENT( magics.unrotate( r ), p );
}

CASE( "test_rotation_angle_only" ) {
    Config config;
    config.set( "rotation_angle", -180. );
    Rotation rotation( config );
    MagicsRotation magics( rotation.southPole() );

    Log::info() << rotation << std::endl;

    EXPECT( rotation.rotated() );

    PointLonLat p, r;

    p = {0., 90.};
    r = {-180., 90};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );

    p = {0., 0.};
    r = {-180., 0.};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );

    p = {270., 25.};
    r = {90., 25.};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );

    p = {-180., 45.};
    r = {-360., 45.};
    EXPECT_EQUIVALENT( rotation.rotate( p ), r );
    EXPECT_EQUIVALENT( rotation.unrotate( r ), p );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
