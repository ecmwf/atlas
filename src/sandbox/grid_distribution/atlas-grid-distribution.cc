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

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
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
    std::string usage() override { return name() + " (--grid.name=name|--grid.json=path) [OPTION]... OUTPUT [--help]"; }

public:
    Tool( int argc, char** argv );

private:
    std::string key;
    PathName path_in;
    PathName path_out;
};

//-----------------------------------------------------------------------------

Tool::Tool( int argc, char** argv ) : AtlasTool( argc, argv ) {
    add_option( new SimpleOption<std::string>(
        "grid.name", "Grid unique identifier\n" + indent() + "     Example values: N80, F40, O24, L32" ) );
    add_option( new SimpleOption<PathName>( "grid.json", "Grid described by json file" ) );
    add_option( new SimpleOption<std::string>( "partitioner", "Partitioner to be used" ) );
    add_option( new SimpleOption<long>( "partitions", "Number of partitions" ) );
}

//-----------------------------------------------------------------------------

int Tool::execute( const Args& args ) {
    key = "";
    args.get( "grid.name", key );

    std::string path_in_str = "";
    if ( args.get( "grid.json", path_in_str ) ) {
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
    else if ( path_in.path().size() ) {
        Log::info() << "Creating grid from file " << path_in << std::endl;
        Log::debug() << Config( path_in ) << std::endl;
        try {
            grid = Grid( Config( path_in ) );
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

    Log::debug() << "Domain: " << grid.domain() << std::endl;
    Log::debug() << "Periodic: " << grid.periodic() << std::endl;

    std::string partitioner_type;
    if ( not args.get( "partitioner", partitioner_type ) ) {
        partitioner_type = "equal_regions";
    }

    long N = mpi::comm().size();
    args.get( "partitions", N );

    if ( mpi::comm().rank() == 0 ) {
        grid::Partitioner partitioner( partitioner_type, N );
        grid::Distribution distribution = partitioner.partition( grid );

        Log::info() << distribution << std::endl;

        std::vector<std::vector<double>> x( N );
        std::vector<std::vector<double>> y( N );
        for ( long p = 0; p < N; ++p ) {
            size_t nb_pts = distribution.nb_pts()[p];
            x[p].reserve( nb_pts );
            y[p].reserve( nb_pts );
        }

        size_t n = 0;
        for ( PointXY pxy : grid.xy() ) {
            size_t p = distribution.partition( n++ );
            x[p].push_back( pxy.x() );
            y[p].push_back( pxy.y() );
        }

        std::ofstream f( "grid-distribution.py", std::ios::trunc );
        f << "\n"
             "import matplotlib.pyplot as plt"
             "\n"
             "from matplotlib.path import Path"
             "\n"
             "import matplotlib.patches as patches"
             "\n"
             ""
             "\n"
             "from itertools import cycle"
             "\n"
             "import matplotlib.cm as cm"
             "\n"
             "import numpy as np"
             "\n"
             "cycol = cycle([cm.Paired(i) for i in "
             "np.linspace(0,1,12,endpoint=True)]).next"
             "\n"
             ""
             "\n"
             "fig = plt.figure()"
             "\n"
             "ax = fig.add_subplot(111,aspect='equal')"
             "\n"
             "";

        for ( long p = 0; p < N; ++p ) {
            f << "\n"
                 "x = [";
            for ( const double& _x : x[p] ) {
                f << _x << ", ";
            }
            f << "]";
            f << "\n"
                 "y = [";
            for ( const double& _y : y[p] ) {
                f << _y << ", ";
            }
            f << "]";
            f << "\n"
                 "c = cycol()";
            f << "\n"
                 "ax.scatter(x, y, color=c, marker='o')";
            f << "\n"
                 "";
        }
        f << "\n"
             "ax.set_xlim(  0-5, 360+5)"
             "\n"
             "ax.set_ylim(-90-5,  90+5)"
             "\n"
             "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
             "\n"
             "ax.set_yticks([-90,-45,0,45,90])"
             "\n"
             "plt.grid()"
             "\n"
             "plt.show()";
    }
    return success();
}

//------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    Tool tool( argc, argv );
    return tool.start();
}
