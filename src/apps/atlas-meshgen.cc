/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildStatistics.h"
#include "atlas/mesh/actions/BuildTorusXYZField.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Bytes.h"
#include "eckit/runtime/Main.h"
#include "eckit/runtime/Tool.h"

//------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::mesh::actions;
using namespace atlas::grid;
using namespace atlas::functionspace;
using namespace atlas::mesh;
using atlas::util::Config;
using eckit::PathName;

//------------------------------------------------------------------------------

class Meshgen2Gmsh : public AtlasTool {
    virtual int execute( const Args& args );
    virtual std::string briefDescription() { return "Mesh generator for Structured compatible meshes"; }
    virtual std::string usage() { return name() + " (--grid.name=name|--grid.json=path) [OPTION]... OUTPUT [--help]"; }

public:
    Meshgen2Gmsh( int argc, char** argv );

private:
    std::string key;
    long halo;
    bool edges;
    bool brick;
    bool stats;
    bool info;
    bool ghost;
    bool binary;
    std::string identifier;
    PathName path_in;
    PathName path_out;
};

//-----------------------------------------------------------------------------

Meshgen2Gmsh::Meshgen2Gmsh( int argc, char** argv ) : AtlasTool( argc, argv ) {
    add_option( new SimpleOption<std::string>(
        "grid.name", "Grid unique identifier\n" + indent() + "     Example values: N80, F40, O24, L32" ) );
    add_option( new SimpleOption<PathName>( "grid.json", "Grid described by json file" ) );
    add_option( new SimpleOption<double>( "angle", "Maximum element-edge slant deviation from meridian in degrees. \n" +
                                                       indent() + "     Value range between 0 and 30\n" + indent() +
                                                       "         0: Mostly triangular, with only perfect quads\n" +
                                                       indent() +
                                                       "        30: Mostly skewed quads with only triags when "
                                                       "skewness becomes too large\n" +
                                                       indent() + "        -1: Only triangles" ) );

    add_option( new SimpleOption<bool>( "include_pole", "Include pole point" ) );
    add_option( new SimpleOption<bool>( "patch_pole", "Patch poles with elements." ) );
    add_option( new SimpleOption<bool>( "ghost", "Output ghost elements" ) );
    add_option( new Separator( "Advanced" ) );
    add_option( new SimpleOption<long>( "halo", "Halo size" ) );
    add_option( new SimpleOption<bool>( "edges", "Build edge datastructure" ) );
    add_option( new SimpleOption<bool>( "brick", "Build brick dual mesh" ) );
    add_option( new SimpleOption<bool>( "stats", "Write statistics file" ) );
    add_option( new SimpleOption<bool>( "info", "Write Info" ) );
    add_option( new SimpleOption<bool>( "binary", "Write binary file" ) );
    add_option( new SimpleOption<std::string>( "generator", "Mesh generator" ) );
    add_option( new SimpleOption<std::string>( "partitioner", "Mesh partitioner" ) );
    add_option( new SimpleOption<bool>( "periodic_x", "periodic mesh in x-direction" ) );
    add_option( new SimpleOption<bool>( "periodic_y", "periodic mesh in y-direction" ) );
    add_option( new SimpleOption<bool>( "torus", "Output mesh as torus" ) );
    add_option( new SimpleOption<bool>( "lonlat", "Output mesh in lon-lat coordinates" ) );
    add_option( new SimpleOption<bool>( "3d",
                                        "Output mesh as sphere, and generate "
                                        "mesh connecting East and West in "
                                        "case serial" ) );
}

//-----------------------------------------------------------------------------

int Meshgen2Gmsh::execute( const Args& args ) {
    key = "";
    args.get( "grid.name", key );

    edges = false;
    args.get( "edges", edges );
    stats = false;
    args.get( "stats", stats );
    info = false;
    args.get( "info", info );
    halo = 0;
    args.get( "halo", halo );
    bool dim_3d = false;
    args.get( "3d", dim_3d );
    brick = false;
    args.get( "brick", brick );
    ghost = false;
    args.get( "ghost", ghost );
    binary = false;
    args.get( "binary", binary );

    std::string path_in_str = "";
    if ( args.get( "grid.json", path_in_str ) ) {
        path_in = path_in_str;
    }

    if ( args.count() ) {
        path_out = args( 0 );
    }
    else {
        path_out = "mesh.msh";
    }

    if ( path_in_str.empty() && key.empty() ) {
        Log::warning() << "missing argument --grid.name or --grid.json" << std::endl;
        Log::warning() << "Usage: " << usage() << std::endl;
        return failed();
    }

    if ( edges ) {
        halo = std::max( halo, 1l );
    }

    StructuredGrid grid;
    if ( key.size() ) {
        try {
            grid = Grid( key );
        }
        catch ( eckit::Exception& ) {
        }
    }
    else if ( path_in.path().size() ) {
        Log::info() << "Creating grid from file " << path_in << std::endl;
        Log::debug() << Config( path_in ) << std::endl;
        try {
            grid = Grid( Config( path_in ) );
        }
        catch ( eckit::Exception& ) {
        }
    }
    else {
        Log::error() << "No grid specified." << std::endl;
    }

    if ( !grid ) {
        return failed();
    }

    Log::debug() << "Domain: " << grid.domain() << std::endl;
    Log::debug() << "Periodic: " << grid.periodic() << std::endl;
    Log::debug() << "Spec: " << grid.spec() << std::endl;

    std::string Implementationype = ( RegularGrid( grid ) ? "regular" : "structured" );
    args.get( "generator", Implementationype );
    eckit::LocalConfiguration meshgenerator_config( args );
    if ( mpi::comm().size() > 1 || edges ) {
        meshgenerator_config.set( "3d", false );
    }

    MeshGenerator meshgenerator( Implementationype, meshgenerator_config );

    Mesh mesh;
    try {
        mesh = meshgenerator.generate( grid );
    }
    catch ( eckit::Exception& e ) {
        Log::error() << e.what() << std::endl;
        Log::error() << e.callStack() << std::endl;
        throw e;
    }

    if ( grid.projection().units() == "degrees" ) {
        functionspace::NodeColumns nodes_fs( mesh, option::halo( halo ) );
    }
    else {
        Log::warning() << "Not yet implemented: building halo's with projections "
                          "not defined in degrees"
                       << std::endl;
        Log::warning() << "units: " << grid.projection().units() << std::endl;
    }
    if ( edges && grid.projection().units() == "degrees" ) {
        functionspace::EdgeColumns edges_fs( mesh, option::halo( halo ) );
        if ( brick ) {
            build_brick_dual_mesh( grid, mesh );
        }
        else {
            build_median_dual_mesh( mesh );
        }
    }

    if ( stats ) {
        build_statistics( mesh );
    }

    bool torus = false;
    args.get( "torus", torus );
    if ( torus ) {
        dim_3d = true;
        Log::debug() << "Building xyz representation for nodes on torus" << std::endl;
        mesh::actions::BuildTorusXYZField( "xyz" )( mesh, grid.domain(), 5., 2., grid.nxmax(), grid.ny() );
    }

    bool lonlat = false;
    args.get( "lonlat", lonlat );

    atlas::output::Gmsh gmsh(
        path_out, Config( "info", info )( "ghost", ghost )( "coordinates", dim_3d ? "xyz" : lonlat ? "lonlat" : "xy" )(
                      "edges", edges )( "binary", binary ) );
    Log::info() << "Writing mesh to gmsh file \"" << path_out << "\" generated from grid \"" << grid.name() << "\""
                << std::endl;
    gmsh.write( mesh );

    if ( info ) {
        Log::info() << "Partitioning graph: \n" << mesh.partitionGraph() << std::endl;
        Log::info() << "Mesh partition footprint: " << eckit::Bytes( mesh.footprint() ) << std::endl;
        for ( idx_t jhalo = 0; jhalo <= halo; ++jhalo ) {
            mesh.polygon( jhalo ).outputPythonScript( "polygon_halo" + std::to_string( jhalo ) + ".py" );
        }
    }
    return success();
}

//------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    Meshgen2Gmsh tool( argc, argv );
    return tool.start();
}
