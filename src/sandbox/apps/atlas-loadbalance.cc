/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>
#include <sstream>

#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

//------------------------------------------------------------------------------------------------------

using eckit::Resource;
using namespace atlas;
using namespace atlas::mesh::actions;
using namespace atlas::grid;
using namespace atlas::functionspace;
using namespace atlas::mesh;

//------------------------------------------------------------------------------------------------------

class AtlasLoadbalance : public eckit::Tool {
    void run() override;

public:
    AtlasLoadbalance( int argc, char** argv ) : eckit::Tool( argc, argv ) {
        bool help = Resource<bool>( "--help", false );

        do_run = true;

        std::string help_str =
            "NAME\n"
            "       atlas-loadbalance - <TODO>\n"
            "\n"
            "SYNOPSIS\n"
            "       atlas-loadbalance GRID [OPTION]... [--help] \n"
            "\n"
            "DESCRIPTION\n"
            "\n"
            "       GRID: unique identifier for grid \n"
            "           Example values: N80, F40, O24, L32\n"
            "\n"
            "       --halo       Output file for mesh\n"
            "\n"
            "AUTHOR\n"
            "       Written by Willem Deconinck.\n"
            "\n"
            "ECMWF                        September 2015";
        if ( help ) {
            atlas::Log::info() << help_str << std::endl;
            do_run = false;
        }

        if ( argc == 1 ) {
            atlas::Log::info() << "usage: atlas-loadbalance GRID [OPTION]... [--help]" << std::endl;
            do_run = false;
        }

        atlas::initialize( argc, argv );

        key = "";
        for ( int i = 0; i < argc; ++i ) {
            if ( i == 1 && argv[i][0] != '-' ) {
                key = std::string( argv[i] );
            }
        }

        halo   = Resource<int>( "--halo", 1 );
        output = Resource<std::string>( "--output", "" );
    }

private:
    bool do_run;
    std::string key;
    int halo;
    std::string output;
    std::string identifier;
};

//------------------------------------------------------------------------------------------------------

void AtlasLoadbalance::run() {
    if ( !do_run ) {
        return;
    }

    StructuredGrid grid;
    try {
        grid = Grid( key );
    }
    catch ( eckit::Exception& err ) {
    }

    if ( !grid ) {
        return;
    }
    MeshGenerator meshgenerator( "structured" );
    Mesh mesh = meshgenerator.generate( grid );

    functionspace::NodeColumns nodes( mesh, option::halo( halo ) );

    if ( output.size() ) {
        write_load_balance_report( mesh, output );
    }
    else {
        std::stringstream s;
        write_load_balance_report( mesh, s );

        if ( mpi::comm().rank() == 0 ) {
            std::cout << s.str() << std::endl;
        }
    }
    atlas::finalize();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    AtlasLoadbalance tool( argc, argv );
    return tool.start();
}
