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

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

//------------------------------------------------------------------------------------------------------

using namespace atlas;
using atlas::util::Config;

//------------------------------------------------------------------------------------------------------

class AtlasLoadbalance : public AtlasTool {
public:
    std::string briefDescription() override { return "Print load balance for meshed grid"; }
    std::string usage() override { return name() + " GRID [OPTION]... [--help]"; }
    std::string longDescription() override {
        return "The 'GRID' argument can be either the name of a named grid, orthe path to a"
               " YAML configuration file that describes the grid.\n"
               "Example values for grid names are: N80, F40, O24, L64x33, CS-ED-12. See the program "
               "'atlas-grids' for a list of named grids.";
    }
    int numberOfPositionalArguments() override { return 1; }
    int minimumPositionalArguments() override { return 1; }

    AtlasLoadbalance(int argc, char** argv): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>(
            "meshgenerator", "Mesh generator [structured,regular,delaunay,cubedsphere] (default depends on GRID)"));
        add_option(new SimpleOption<std::string>(
            "partitioner",
            "Grid partitioner [equal_regions,checkerboard,equal_bands,regular_bands] (default depends on GRID)"));
        add_option(new Separator("Advanced"));
        add_option(new SimpleOption<long>("halo", "Halo size"));
    }

    int execute(const Args& args) override {
        std::string key = get_positional_arg(args, 0);

        Grid grid;
        if (key.size()) {
            eckit::PathName path_in{key};
            try {
                if (path_in.exists()) {
                    Log::info() << "Creating grid from file " << path_in << std::endl;
                    Log::debug() << Config(path_in) << std::endl;
                    grid = Grid(Config{path_in});
                }
                else {
                    Log::info() << "Creating grid from name " << key << std::endl;
                    grid = Grid(key);
                }
            }
            catch (eckit::Exception& e) {
                Log::error() << e.what() << std::endl;
                if (path_in.exists()) {
                    Log::error() << "Could not generate mesh for grid defined in file \"" << path_in << "\""
                                 << std::endl;
                }
                else {
                    Log::error() << "Could not generate mesh for grid \"" << key << "\"" << std::endl;
                }
                return failed();
            }
        }
        else {
            Log::error() << "No grid specified." << std::endl;
            return failed();
        }


        auto meshgenerator = make_meshgenerator(grid, args);
        auto partitioner   = make_partitioner(grid, args);
        int halo           = args.getInt("halo", 0);

        Log::info() << "Input:" << std::endl;
        Log::info() << "- grid:          " << grid.name() << std::endl;
        Log::info() << "- meshgenerator: " << meshgenerator.type() << std::endl;
        Log::info() << "- partitioner:   " << partitioner.type() << std::endl;
        Log::info() << "- halo:          " << halo << std::endl;


        auto mesh = meshgenerator.generate(grid, partitioner);
        functionspace::NodeColumns nodes(mesh, option::halo(halo));

        {
            std::stringstream s;
            mesh::actions::write_load_balance_report(mesh, s);

            if (mpi::comm().rank() == 0) {
                std::cout << s.str() << std::endl;
            }
        }
        return success();
    }

    MeshGenerator make_meshgenerator(const Grid& grid, const Args& args) {
        auto config = grid.meshgenerator();  // recommended by the grid itself
        if (args.has("meshgenerator")) {
            config.set("type", args.getString("meshgenerator"));
        }
        if (args.has("halo")) {
            config.set("halo", args.getInt("halo"));
        }
        return MeshGenerator{config};
    }

    grid::Partitioner make_partitioner(const Grid& grid, const Args& args) {
        auto config = grid.partitioner();  // recommended by the grid itself
        if (args.has("partitioner")) {
            config.set("type", args.getString("partitioner"));
        }
        if (args.has("regular")) {
            config.set("regular", args.getBool("regular"));
        }
        if (args.has("partitions")) {
            config.set("partitions", args.getInt("partitions"));
        }
        return grid::Partitioner{config};
    }
};

//------------------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    AtlasLoadbalance tool(argc, argv);
    return tool.start();
}
