/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstdlib>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::PointCloud;
using atlas::functionspace::StructuredColumns;

namespace atlas {
namespace test {


util::Config scheme() {
    util::Config scheme;
    scheme.set( "type", "structured-linear2D" );
    scheme.set( "experimental-mpi-interpolation", true );
    return scheme;
}

std::string input_gridname( const std::string& default_grid ) {
    return eckit::Resource<std::string>( "--input-grid", default_grid );
}

double vortex_rollup( double lon, double lat, double t ) {
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
    if ( rho != 0. ) {
        omega = 0.5 * 3 * std::sqrt( 3 ) * a * Omega * sqr( sech( rho ) ) * std::tanh( rho ) / rho;
    }
    double q = 1. - std::tanh( 0.2 * rho * std::sin( lambda_prime - omega / a * t ) );
    return q;
};


FunctionSpace output_functionspace_match() {
    std::vector<PointXY> points;
    if ( mpi::size() == 2 ) {
        if ( mpi::rank() == 0 ) {
            points = std::vector<PointXY>{
                {45., 45.}, {90., 45.}, {135., 45.}, {180., 45.}, {225., 45.}, {270., 45.}, {315., 45.},
            };
        }
        if ( mpi::rank() == 1 ) {
            points = std::vector<PointXY>{
                {45., -45.}, {90., -45.}, {135., -45.}, {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
            };
        }
    }
    else if ( mpi::size() == 1 ) {
        points = std::vector<PointXY>{
            {45., 45.},  {90., 45.},  {135., 45.},  {180., 45.},  {225., 45.},  {270., 45.},  {315., 45.},
            {45., -45.}, {90., -45.}, {135., -45.}, {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
        };
    }
    else {
        return FunctionSpace();
    }
    return PointCloud( points );
}


FunctionSpace output_functionspace_nomatch() {
    std::vector<PointXY> points;
    if ( mpi::size() == 2 ) {
        if ( mpi::rank() == 0 ) {
            points = std::vector<PointXY>{
                {45., 45.}, {90., 45.}, {135., 45.}, {45., -45.}, {90., -45.}, {135., -45.},
            };
        }
        if ( mpi::rank() == 1 ) {
            points = std::vector<PointXY>{
                {180., 45.},  {225., 45.},  {270., 45.},  {315., 45.},
                {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
            };
        }
        else if ( mpi::size() == 1 ) {
            points = std::vector<PointXY>{
                {45., 45.},  {90., 45.},  {135., 45.},  {180., 45.},  {225., 45.},  {270., 45.},  {315., 45.},
                {45., -45.}, {90., -45.}, {135., -45.}, {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
            };
        }
    }
    else {
        return FunctionSpace();
    }
    return PointCloud( points );
}


FunctionSpace output_functionspace_same_points( const StructuredColumns& fs ) {
    std::vector<PointXY> points;
    points.reserve( fs.size() );

    auto lonlat = array::make_view<double, 2>( fs.lonlat() );
    auto ghost  = array::make_view<int, 1>( fs.ghost() );
    for ( idx_t i = 0; i < fs.size(); ++i ) {
        if ( !ghost( i ) ) {
            points.emplace_back( lonlat( i, LON ), lonlat( i, LAT ) );
        }
    }

    return PointCloud( points );
}


FieldSet create_source_fields_vortex( const StructuredColumns& fs, idx_t nb_fields, idx_t nb_levels ) {
    using Value = double;

    FieldSet fields_source;
    auto lonlat = array::make_view<double, 2>( fs.xy() );

    for ( idx_t f = 0; f < nb_fields; ++f ) {
        auto field_source = fields_source.add( fs.createField<Value>() );
        auto source       = array::make_view<Value, 2>( field_source );
        for ( idx_t n = 0; n < fs.size(); ++n ) {
            for ( idx_t k = 0; k < nb_levels; ++k ) {
                source( n, k ) = vortex_rollup( lonlat( n, LON ), lonlat( n, LAT ), 0.5 + double( k ) / 2 );
            }
        }
    }

    return fields_source;
}


FieldSet create_source_fields_rand( const StructuredColumns& fs, idx_t nb_fields, idx_t nb_levels ) {
    using Value = double;

    FieldSet fields_source;

    for ( idx_t f = 0; f < nb_fields; ++f ) {
        auto field  = fields_source.add( fs.createField<Value>( option::levels( nb_levels ) ) );
        auto source = array::make_view<Value, 2>( field );
        for ( idx_t n = 0; n < field.shape( 0 ); ++n ) {
            for ( idx_t l = 0; l < nb_levels; ++l ) {
                source( n, l ) = Value( std::rand() );
            }
        }
    }

    return fields_source;
}


FieldSet create_source_fields_lonlat( const StructuredColumns& fs ) {
    using Value = double;

    FieldSet fields_source;
    auto lonlat = array::make_view<double, 2>( fs.xy() );

    for ( auto ll : {LON, LAT} ) {
        auto field  = fields_source.add( fs.createField<Value>() );
        auto source = array::make_view<Value, 2>( field );
        for ( idx_t n = 0; n < fs.xy().shape( 0 ); ++n ) {
            source( n, 0 ) = lonlat( n, ll );
        }
    }

    return fields_source;
}


FieldSet create_target_fields( FunctionSpace& fs, idx_t nb_fields, idx_t nb_levels ) {
    using Value = double;
    FieldSet fields_target;
    for ( idx_t f = 0; f < nb_fields; ++f ) {
        fields_target.add( fs.createField<Value>( option::levels( nb_levels ) ) );
    }
    return fields_target;
}


CASE( "test_match" ) {
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid( input_gridname( "O32" ) );
    StructuredColumns input_fs( input_grid, option::halo( 1 ) | option::levels( nb_levels ) );

    FunctionSpace output_fs = output_functionspace_match();

    Interpolation interpolation( scheme(), input_fs, output_fs );

    FieldSet fields_source = create_source_fields_vortex( input_fs, nb_fields, nb_levels );
    FieldSet fields_target = create_target_fields( output_fs, nb_fields, nb_levels );

    interpolation.execute( fields_source, fields_target );
}

CASE( "test_nomatch" ) {
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid( input_gridname( "O32" ) );
    StructuredColumns input_fs( input_grid, option::halo( 1 ) | option::levels( nb_levels ) );

    FunctionSpace output_fs = output_functionspace_nomatch();

    FieldSet fields_source = create_source_fields_vortex( input_fs, nb_fields, nb_levels );
    FieldSet fields_target = create_target_fields( output_fs, nb_fields, nb_levels );

    Interpolation interpolation( scheme(), input_fs, output_fs );
    interpolation.execute( fields_source, fields_target );
}


CASE( "test_nomatch_same_points" ) {
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid( input_gridname( "O32" ) );

    StructuredColumns input_fs( input_grid, option::halo( 1 ) );
    FieldSet fields_source = create_source_fields_rand( input_fs, nb_fields, nb_levels );

    FunctionSpace output_fs = output_functionspace_same_points( input_fs );
    FieldSet fields_target  = create_target_fields( output_fs, nb_fields, nb_levels );

    ASSERT( !fields_source.empty() && !fields_target.empty() );
    ASSERT( input_fs.size() >= output_fs.size() );

    Interpolation interpolation( scheme(), input_fs, output_fs );
    interpolation.execute( fields_source, fields_target );

    double eps = 1e-3;

    ASSERT( fields_source.size() == fields_target.size() );
    for ( auto f = 0; f < fields_target.size(); ++f ) {
        const auto source_values = array::make_view<double, 2>( fields_source[f] );
        const auto target_values = array::make_view<double, 2>( fields_target[f] );

        ASSERT( source_values.shape( 0 ) >= target_values.shape( 0 ) );
        ASSERT( source_values.shape( 1 ) == target_values.shape( 1 ) );
        auto ghost = array::make_view<int, 1>( input_fs.ghost() );

        for ( idx_t i = 0; i < source_values.shape( 0 ); ++i ) {
            if ( !ghost( i ) ) {
                for ( idx_t l = 0; l < source_values.shape( 1 ); ++l ) {
                    EXPECT( is_approximately_equal( source_values( i, l ), target_values( i, l ), eps ) );
                }
            }
        }
    }
}


CASE( "test_nomatch_lonlat" ) {
    Grid input_grid( input_gridname( "O90" ) );  // 1-degree/1-degree grid approx.

    StructuredColumns input_fs( input_grid, option::halo( 1 ) | option::levels( 1 ) );
    FieldSet fields_source = create_source_fields_lonlat( input_fs );
    ASSERT( fields_source.size() == 2 );

    FunctionSpace output_fs = output_functionspace_nomatch();
    FieldSet fields_target  = create_target_fields( output_fs, fields_source.size(), 1 );
    ASSERT( fields_source.size() == fields_target.size() );

    Interpolation interpolation( scheme(), input_fs, output_fs );
    interpolation.execute( fields_source, fields_target );

    const auto lon = array::make_view<double, 2>( fields_target[LON] );
    const auto lat = array::make_view<double, 2>( fields_target[LAT] );
    ASSERT( lon.size() == lat.size() );

    double eps  = 1e-3;
    auto lonlat = array::make_view<double, 2>( output_fs.lonlat() );
    for ( auto i = 0; i < lon.size(); ++i ) {
        EXPECT( is_approximately_equal( lonlat( i, LON ), lon( i, 0 ), eps ) );
        EXPECT( is_approximately_equal( lonlat( i, LAT ), lat( i, 0 ), eps ) );
    }
}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
