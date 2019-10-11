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
#include <string>
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "eckit/config/Resource.h"
#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/log/Plural.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatPolygon.h"


using namespace atlas;


class AtlasParallelInterpolation : public AtlasTool {
    int execute( const AtlasTool::Args& args ) override;
    std::string briefDescription() override { return "Demonstration of parallel interpolation"; }
    std::string usage() override {
        return name() +
               " [--source-gridname=gridname] "
               "[--target-gridname=gridname] [OPTION]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasParallelInterpolation( int argc, char* argv[] ) : AtlasTool( argc, argv ) {

        add_option( new SimpleOption<std::string>( "source-gridname", "source gridname" ) );
        add_option( new SimpleOption<std::string>( "target-gridname", "target gridname" ) );

        add_option( new SimpleOption<std::string>( "method", "interpolation method (default linear)" ) );
        add_option( new SimpleOption<bool>( "matrix-free", "Don't store matrix for consecutive interpolations" ) );

        add_option( new SimpleOption<std::string>( "partitioner",
                                                   "source partitioner (equal_regions (default), ...)" ) );
        
        add_option( new SimpleOption<bool>( "output-gmsh",
                                            "Output gmsh files src_field.msh and tgt_field.msh" ) );
        add_option(
            new SimpleOption<bool>( "output-polygons", "Output Python script that plots partitions polygons" ) );
        
        add_option( new SimpleOption<long>( "vortex-rollup","Value that controls vortex rollup (default = 0)") );
        add_option( new SimpleOption<bool>( "with-backwards","Do backwards interpolation") );
        
    }
};

static Config processed_config( const eckit::Configuration& _config ) {
    Config config;
    if( _config.has("partitioner") ) {
        config.set("partitioner",option::type(_config.getString("partitioner")));
    }
    std::string scheme_str = _config.getString("method","linear");
    if ( scheme_str == "linear" ) {
        config.set( "type", "structured-linear2D" );
        config.set( "halo", 1 );
        // The stencil does not require any halo, but we set it to 1 for pole treatment!
    }
    if ( scheme_str == "cubic" ) {
        config.set( "type", "structured-cubic2D" );
        config.set( "halo", 2 );
    }
    if ( scheme_str == "quasicubic" ) {
        config.set( "type", "structured-quasicubic2D" );
        config.set( "halo", 2 );
    }
    config.set( "name", scheme_str );
    config.set("matrix_free",_config.getBool("matrix-free",false));
    return config;
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

grid::Distribution distribution( Grid& grid, FunctionSpace& src ) {
    ATLAS_TRACE("Computing distribution from source");
    functionspace::detail::PartitionPolygon p(*src.get(),0);
    
    int rank = mpi::comm().rank();
    util::LonLatPolygon poly{ p.lonlat() };
    std::vector<int> part; part.reserve( grid.size() );
    for( auto p : grid.lonlat() ) {
        if( poly.contains(p) ) {
            part.emplace_back( rank );
        } 
        else {
            part.emplace_back( -1 );
        }
    }
    mpi::comm().allReduceInPlace(part.begin(),part.end(),eckit::mpi::max());
    return grid::Distribution(part.size(),part.data());
}


int AtlasParallelInterpolation::execute( const AtlasTool::Args& args ) {

    ATLAS_TRACE("AtlasParallelInterpolation::execute");
    auto source_gridname = args.getString( "source-gridname", "O32" );
    auto target_gridname = args.getString( "target-gridname", "O64" );
    auto config = processed_config(args);

    Log::info() << "atlas-parallel-interpolation from source grid " << source_gridname << " to " << target_gridname
                << std::endl;
    Grid src_grid( source_gridname );
    Grid tgt_grid( target_gridname );

    int nlev = 0;

    auto src_fs  = functionspace::StructuredColumns{ src_grid, config | option::levels( nlev ) };

    if( args.getBool("output-polygons",false) ) {
        functionspace::detail::PartitionPolygon src_poly(*src_fs.get(),0);
        src_poly.outputPythonScript("src-polygons.py");
    }
    
    
    auto tgt_dist = distribution(tgt_grid,src_fs);
    auto tgt_fs = functionspace::StructuredColumns{ tgt_grid, tgt_dist, config | option::levels( nlev ) };

    if( args.getBool("output-polygons",false) ) {
        functionspace::detail::PartitionPolygon tgt_poly(*tgt_fs.get(),0);
        tgt_poly.outputPythonScript("tgt-polygons.py",Config("nodes",false));
    }

    
    Field src_field = src_fs.createField<double>( option::name( "source" ) );
    Field tgt_field = tgt_fs.createField<double>( option::name( "target" )  );

    // Initialize source
    {
        ATLAS_TRACE("Initialize source");
        auto lonlat = array::make_view<double, 2>( src_fs.xy() );
        auto source = array::make_view<double, 1>( src_field );
        idx_t size = src_fs.size();
        idx_t k = args.getInt("vortex-rollup",0);
        for ( idx_t n = 0; n < size; ++n ) {
            source( n ) = vortex_rollup( lonlat( n, LON ), lonlat( n, LAT ), 0.5 + double( k ) / 2 );
        };
    }
    

    src_field.haloExchange();
    
    output::Gmsh src_gmsh("src_field.msh");
    output::Gmsh tgt_gmsh("tgt_field.msh",Config("ghost",true));
    
    if( args.getBool("output-gmsh",false) ) {
        src_gmsh.write( src_field );
    }
    
    {
        ATLAS_TRACE("Interpolation: Source to Target");
        Interpolation interpolation_fwd( config, src_fs, tgt_fs );
        interpolation_fwd.execute( src_field, tgt_field );
    }    

    tgt_field.haloExchange();
    if( args.getBool("output-gmsh",false) ) {
        tgt_gmsh.write( tgt_field );
    }
    
    if( args.getBool("with-backwards",false))
    {
       {
            ATLAS_TRACE("Interpolation: Target to Source");
            Interpolation interpolation_bwd( config, tgt_fs, src_fs );
            interpolation_bwd.execute( tgt_field, src_field );
        }
        
        if( args.getBool("output-gmsh",false) ) {
            src_gmsh.write( src_field );
        }
    }
    
    return success();
}


int main( int argc, char* argv[] ) {
    AtlasParallelInterpolation tool( argc, argv );
    return tool.start();
}
