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
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/meshgenerator.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

//------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;
using atlas::util::Config;
using eckit::PathName;

//------------------------------------------------------------------------------

class Tool : public AtlasTool {
    int execute( const Args& args ) override;
    std::string briefDescription() override {
        return "Tool to generate a python script that plots the grid-distribution "
               "of a given grid";
    }
    std::string usage() override { return name() + " --grid=name [OPTION]... OUTPUT [--help]"; }

public:
    Tool( int argc, char** argv );

private:
    std::string key;
    PathName path_in;
    PathName path_out;
};

//-----------------------------------------------------------------------------

static int halo_default() {
    return 2;
}

Tool::Tool( int argc, char** argv ) : AtlasTool( argc, argv ) {
    add_option( new SimpleOption<std::string>(
        "grid", "Grid unique identifier\n" + indent() + "     Example values: N80, F40, O24, L32" ) );
    add_option(
        new SimpleOption<long>( "halo", "Number of halos (default=" + std::to_string( halo_default() ) + ")" ) );
}

//-----------------------------------------------------------------------------

int Tool::execute( const Args& args ) {
    Trace timer( Here(), displayName() );
    key = "";
    args.get( "grid", key );

    std::string path_in_str = "";
    if ( args.get( "grid", path_in_str ) ) {
        path_in = path_in_str;
    }

    StructuredGrid grid;
    if ( key.size() ) {
        try {
            grid = Grid( key );
        }
        catch ( eckit::Exception& e ) {
        }
    }
    else {
        Log::error() << "No grid specified." << std::endl;
    }

    if ( !grid ) {
        return failed();
    }

    size_t halo = args.getLong( "halo", halo_default() );

    Log::info() << "Configuration" << std::endl;
    Log::info() << "~~~~~~~~~~~~~" << std::endl;
    Log::info() << "  Grid   : " << grid.name() << std::endl;
    Log::info() << "  Halo   : " << halo << std::endl;
    Log::info() << "  MPI    : " << mpi::comm().size() << std::endl;
    Log::info() << "  OpenMP : " << atlas_omp_get_max_threads() << std::endl;

    MeshGenerator meshgenerator( "structured", util::Config( "partitioner", "equal_regions" ) );

    size_t iterations = 1;
    for ( size_t i = 0; i < iterations; ++i ) {
        ATLAS_TRACE( "iteration" );
        Mesh mesh = meshgenerator.generate( grid );
        numerics::fvm::Method fvm( mesh, mesh::Halo( halo ) );
        // mesh::actions::build_halo( mesh, halo );
        mpi::comm().barrier();
    }
    timer.stop();
    Log::info() << Trace::report() << std::endl;
    return success();
}

//------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    Tool tool( argc, argv );
    return tool.start();
}
