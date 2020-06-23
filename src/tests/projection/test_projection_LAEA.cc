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
#include <vector>
#include "atlas/grid.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"


using atlas::util::Config;
using ListPointXY = std::vector<atlas::PointXY>;
namespace atlas {
namespace test {

static const double km = 1000.;

//-----------------------------------------------------------------------------

CASE( "test_projection_LAEA" ) {
    auto points_xy  = ListPointXY{// in Meters!!!
                                 {0., 0.},
                                 {100. * km, 0.},
                                 {0., 100. * km}};
    auto projection = Projection{
        Config( "type", "lambert_azimuthal_equal_area" )( "central_longitude", -67. )( "standard_parallel", 50. )};
    EXPECT( projection.type() == "lambert_azimuthal_equal_area" );
    EXPECT( projection.units() == "meters" );
    EXPECT( projection.strictlyRegional() );

    for ( auto pxy : points_xy ) {
        double tolerance_micrometers = 1.e-6;
        auto p_identity              = projection.xy( projection.lonlat( pxy ) );
        EXPECT( is_approximately_equal( p_identity.x(), pxy.x(), tolerance_micrometers ) );
        EXPECT( is_approximately_equal( p_identity.y(), pxy.y(), tolerance_micrometers ) );
    }
}

//-----------------------------------------------------------------------------

CASE( "test_grid_creation_from_GRIB" ) {
    auto mir_parametrisation = []() -> Config {
        Config parametrisation;
        parametrisation.set( "gridType", "lambert_azimuthal_equal_area" );
        parametrisation.set( "standardParallelInDegrees", 50. );
        parametrisation.set( "centralLongitudeInDegrees", -67. );
        parametrisation.set( "latitudeOfFirstGridPointInDegrees", 75. );
        parametrisation.set( "longitudeOfFirstGridPointInDegrees", -90. );
        parametrisation.set( "xDirectionGridLengthInMetres", 100. * km );
        parametrisation.set( "yDirectionGridLengthInMetres", 10. * km );
        parametrisation.set( "numberOfPointsAlongXAxis", 100 );
        parametrisation.set( "numberOfPointsAlongYAxis", 50 );
        parametrisation.set( "jScansPositively", 0L );
        parametrisation.set( "radius", util::Earth::radius() );
        return parametrisation;
    };

    auto from_GRIB = []( const eckit::Parametrisation& parametrisation ) -> Grid::Spec {
        // Convert parametrisation to Grid::Spec
        std::vector<double> firstPointLonLat( 2 );
        double standard_parallel, central_longitude, dx, dy;
        long nx, ny;
        double radius         = util::Earth::radius();
        long jScansPositively = 0;

        parametrisation.get( "standardParallelInDegrees", standard_parallel );
        parametrisation.get( "centralLongitudeInDegrees", central_longitude );
        parametrisation.get( "radius", radius );
        parametrisation.get( "longitudeOfFirstGridPointInDegrees", firstPointLonLat[0] );
        parametrisation.get( "latitudeOfFirstGridPointInDegrees", firstPointLonLat[1] );
        parametrisation.get( "numberOfPointsAlongXAxis", nx );
        parametrisation.get( "numberOfPointsAlongYAxis", ny );
        parametrisation.get( "xDirectionGridLengthInMetres", dx );
        parametrisation.get( "yDirectionGridLengthInMetres", dy );
        parametrisation.get( "jScansPositively", jScansPositively );
        std::string firstPointLoc = jScansPositively ? "lonlat(xmin,ymin)" : "lonlat(xmin,ymax)";

        Grid::Spec gridspec;
        gridspec.set( "type", "regional" );  // --> indicates following values will be parsed
        gridspec.set( "nx", nx );
        gridspec.set( "ny", ny );
        gridspec.set( "dx", dx );
        gridspec.set( "dy", dy );
        gridspec.set( firstPointLoc, firstPointLonLat );
        gridspec.set( "y_numbering", -1 );  // Always decrease y for MIR
        gridspec.set( "projection", [&]() {
            Grid::Spec projection;
            projection.set( "type", "lambert_azimuthal_equal_area" );
            projection.set( "standard_parallel", standard_parallel );
            projection.set( "central_longitude", central_longitude );
            projection.set( "radius", radius );
            return projection;
        }() );
        return gridspec;
    };

    // Create the grid from mir_parametrisation
    auto g = Grid( from_GRIB( mir_parametrisation() ) );

    // Compute first point and compare with mir_parametrisation
    {
        auto domain       = RectangularDomain( g.domain() );
        auto firstPointXY = mir_parametrisation().getLong( "jScansPositively", 0L )
                                ? PointXY{domain.xmin(), domain.ymin()}
                                : PointXY{domain.xmin(), domain.ymax()};
        auto firstPointLonLat = g.projection().lonlat( firstPointXY );
        EXPECT( is_approximately_equal( firstPointLonLat.lon(),
                                        mir_parametrisation().getDouble( "longitudeOfFirstGridPointInDegrees" ) ) );
        EXPECT( is_approximately_equal( firstPointLonLat.lat(),
                                        mir_parametrisation().getDouble( "latitudeOfFirstGridPointInDegrees" ) ) );
    }

    // Check if iteration is correct
    {
        auto domain           = RectangularDomain( g.domain() );
        auto firstPointLonLat = g.projection().lonlat( {domain.xmin(), domain.ymax()} );
        auto lastPointLonLat  = g.projection().lonlat( {domain.xmax(), domain.ymin()} );
        long n                = 0;
        for ( PointLonLat p : g.lonlat() ) {
            if ( n == 0 ) {
                EXPECT( is_approximately_equal( p.lon(), firstPointLonLat.lon() ) );
                EXPECT( is_approximately_equal( p.lat(), firstPointLonLat.lat() ) );
            }
            if ( n == g.size() - 1 ) {
                EXPECT( is_approximately_equal( p.lon(), lastPointLonLat.lon() ) );
                EXPECT( is_approximately_equal( p.lat(), lastPointLonLat.lat() ) );
            }
            ++n;
        }
        EXPECT( n == g.size() );
    }

    // Check if bounding box is correct
    {
        RectangularLonLatDomain bb{g.lonlatBoundingBox()};
        const double tolerance = 1.e-6;
        EXPECT_APPROX_EQ( bb.west(), -90., tolerance );
        EXPECT_APPROX_EQ( bb.east(), 41.882046, tolerance );
        EXPECT_APPROX_EQ( bb.south(), 3.818356, tolerance );
        EXPECT_APPROX_EQ( bb.north(), 76.040387, tolerance );
        for ( PointLonLat p : g.lonlat() ) {
            EXPECT( bb.contains( p ) );
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
