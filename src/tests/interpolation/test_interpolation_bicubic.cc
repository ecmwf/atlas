/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/interpolation.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::string input_gridname( const std::string& default_grid ) {
    return eckit::Resource<std::string>( "--input-grid", default_grid );
}

std::string output_gridname( const std::string& default_grid ) {
    return eckit::Resource<std::string>( "--output-grid", default_grid );
}

Grid rotated_mercator() {
    util::Config gridspec;
    gridspec.set( "type", "regional" );
    gridspec.set( "nx", 50 );
    gridspec.set( "ny", 40 );
    gridspec.set( "dx", 50000 );
    gridspec.set( "dy", 50000 );
    gridspec.set( "lonlat(centre)", std::vector<double>{4., 50} );
    gridspec.set( "projection", [] {
        util::Config projection;
        projection.set( "type", "rotated_mercator" );
        projection.set( "north_pole", std::vector<double>{-176., 40.} );
        return projection;
    }() );
    return Grid{gridspec};
}

Grid lambert() {
    util::Config gridspec;
    gridspec.set( "type", "regional" );
    gridspec.set( "nx", 50 );
    gridspec.set( "ny", 40 );
    gridspec.set( "dx", 50000 );
    gridspec.set( "dy", 50000 );
    gridspec.set( "lonlat(centre)", std::vector<double>{4., 50} );
    gridspec.set( "projection", [] {
        util::Config projection;
        projection.set( "type", "lambert" );
        projection.set( "latitude1", 50. );
        projection.set( "longitude0", 4. );
        return projection;
    }() );
    return Grid{gridspec};
}

Grid rotated( const std::string& name ) {
    util::Config gridspec;
    gridspec.set( "name", name );
    gridspec.set( "projection", [] {
        util::Config projection;
        projection.set( "type", "rotated_lonlat" );
        projection.set( "north_pole", std::vector<double>{-176., 40.} );
        return projection;
    }() );
    return Grid{gridspec};
}

FunctionSpace output_functionspace( const Grid& grid ) {
    MeshGenerator meshgen( "structured" );
    Mesh output_mesh = meshgen.generate( grid );
    return NodeColumns{output_mesh};
}

auto vortex_rollup = []( double lon, double lat, double t ) {
    // lon and lat in degrees!

    // Formula found in "A Lagrangian Particle Method with Remeshing for Tracer Transport on the Sphere"
    // by Peter Bosler, James Kent, Robert Krasny, CHristiane Jablonowski, JCP 2015

    lon *= M_PI / 180.;
    lat *= M_PI / 180.;

    auto sqr           = []( const double x ) { return x * x; };
    auto sech          = []( const double x ) { return 1. / std::cosh( x ); };
    const double T     = 1.;
    const double Omega = 2. * M_PI / T;
    t *= T;
    const double lambda_prime = std::atan2( -std::cos( lon - Omega * t ), std::tan( lat ) );
    const double rho          = 3. * std::sqrt( 1. - sqr( std::cos( lat ) ) * sqr( std::sin( lon - Omega * t ) ) );
    double omega              = 0.;
    double a                  = util::Earth::radius();
    if ( rho != 0. ) { omega = 0.5 * 3 * std::sqrt( 3 ) * a * Omega * sqr( sech( rho ) ) * std::tanh( rho ) / rho; }
    double q = 1. - std::tanh( 0.2 * rho * std::sin( lambda_prime - omega / a * t ) );
    return q;
};

CASE( "test_interpolation_cubic_structured" ) {
    Grid grid( input_gridname( "O32" ) );
    StructuredColumns fs( grid, option::halo( 2 ) );

    auto test = [&]( const FunctionSpace& output_fs ) {

        Interpolation interpolation( option::type( "structured-bicubic" ), fs, output_fs );

        Field field_source = fs.createField<double>( option::name( "source" ) );
        Field field_target = output_fs.createField<double>( option::name( "target" ) );

        auto lonlat = array::make_view<double, 2>( fs.xy() );
        auto source = array::make_view<double, 1>( field_source );
        for ( idx_t n = 0; n < fs.size(); ++n ) {
            source( n ) = vortex_rollup( lonlat( n, LON ), lonlat( n, LAT ), 1. );
        }

        interpolation.execute( field_source, field_target );

        output_fs.haloExchange( field_target );

        output::Gmsh gmsh( "cubic-output-section" + std::to_string( _subsection ) + ".msh",
                           util::Config( "coordinates", "xy" ) );
        gmsh.write( NodeColumns( output_fs ).mesh() );
        gmsh.write( field_target );

    };

    SECTION( "Interpolate from " + grid.name() + " to " + output_gridname( "O64" ) ) {
        EXPECT_NO_THROW( test( output_functionspace( Grid{output_gridname( "O64" )} ) ) );
    }
    SECTION( "Interpolate from " + grid.name() + " to rotated " + output_gridname( "O64" ) ) {
        EXPECT_NO_THROW( test( output_functionspace( rotated( output_gridname( "O64" ) ) ) ) );
    }
    SECTION( "Interpolate from " + grid.name() + " to lambert" ) {
        EXPECT_NO_THROW( test( output_functionspace( lambert() ) ) );
    }
    SECTION( "Interpolate from " + grid.name() + " to rotaded_mercator" ) {
        EXPECT_NO_THROW( test( output_functionspace( rotated_mercator() ) ) );
    }
}  // namespace test


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
