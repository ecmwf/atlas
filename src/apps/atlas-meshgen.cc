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
#include "atlas/library.h"
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

MeshGenerator make_meshgenerator(const Grid& grid, const AtlasTool::Args& args) {
    auto config = grid.meshgenerator();  // recommended by the grid itself
    if (args.has("generator")) {
        config.set("type", args.getString("generator"));
    }

    if (args.has("3d")) {
        config.set("3d", true);
    }

    if (mpi::comm().size() > 1 || args.getBool("edges", false)) {
        config.set("3d", false);
    }

    config.set("include_pole", args.getBool("include-pole", false));
    config.set("patch_pole", args.getBool("patch-pole", false));
    if (args.has("fixup")) {
        config.set("fixup", args.getBool("fixup"));
    }

    if (args.has("partition")) {
        config.set("partition", args.getInt("partition"));
        config.set("part", args.getInt("partition"));
    }
    if (args.has("partitions")) {
        config.set("partitions", args.getInt("partitions"));
        config.set("nb_parts", args.getInt("partitions"));
    }
    if (args.has("halo")) {
        config.set("halo", args.getInt("halo"));
    }


    return MeshGenerator{config};
}

Partitioner make_partitioner(const Grid& grid, const AtlasTool::Args& args) {
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
    return Partitioner{config};
}

//------------------------------------------------------------------------------

class Meshgen2Gmsh : public AtlasTool {
    int execute(const Args& args) override;
    std::string briefDescription() override { return "Mesh generator and output to Gmsh format"; }
    std::string usage() override { return name() + " GRID [OUTPUT] [OPTION]... [--help]"; }
    std::string longDescription() override {
        return "    The 'GRID' argument can be either the name of a named grid, orthe path to a"
               " YAML configuration file that describes the grid.\n"
               "Example values for grid names are: N80, F40, O24, L64x33, CS-ED-12. See the program "
               "'atlas-grids' for a list of named grids.\n"
               "\n"
               "    The optional 'OUTPUT' argument contains the path to the output file. "
               "If not given, a default value 'mesh.msh' is employed.";
    }

public:
    Meshgen2Gmsh(int argc, char** argv);

private:
    std::string key;
    long halo;
    bool nodes;
    bool edges;
    bool cells;
    bool brick;
    bool stats;
    bool info;
    bool ghost;
    std::string identifier;
    PathName path_in;
    PathName path_out;
};

//-----------------------------------------------------------------------------

Meshgen2Gmsh::Meshgen2Gmsh(int argc, char** argv): AtlasTool(argc, argv) {
    add_option(new SimpleOption<std::string>("coordinates", "Output mesh in given coordinates"));
    add_option(
        new SimpleOption<bool>("lonlat", "Output mesh in lon,lat coordinates (shorthand for --coordinates=lonlat)"));
    add_option(new SimpleOption<bool>("ij", "Output mesh in i,j coordinates (shorthand for --coordinates=ij)"));
    add_option(new SimpleOption<bool>("3d",
                                      "Output mesh as sphere, and generate "
                                      "mesh connecting East and West in "
                                      "case serial"));
    add_option(new SimpleOption<bool>("ghost", "Output ghost elements"));
    add_option(new SimpleOption<std::string>(
        "generator", "Mesh generator [structured,regular,delaunay,cubedsphere] (default = structured)"));
    add_option(new SimpleOption<std::string>("partitioner",
                                             "Mesh partitioner [equal_regions,checkerboard,equal_bands,regular_bands"));

    add_option(new Separator("Options for `--generator=structured`"));
    add_option(new SimpleOption<bool>("include-pole", "Include pole point"));
    add_option(new SimpleOption<bool>("patch-pole", "Patch poles with elements."));
    add_option(new SimpleOption<double>(
        "angle", "Maximum element-edge slant deviation from meridian in degrees. \n" + indent() +
                     "     Value range between 0 and 30\n" + indent() +
                     "         0: Mostly triangular, with only perfect quads (=default)\n" + indent() +
                     "        30: Mostly skewed quads with only triags when skewness becomes too large\n" + indent() +
                     "        -1: Only triangles"));

    add_option(new Separator("Options for `--partitioner=checkerboard`"));
    add_option(new SimpleOption<bool>("regular", "regular checkerboard partitioner"));

    add_option(new Separator("Advanced"));
    add_option(new SimpleOption<long>("halo", "Halo size"));
    add_option(new SimpleOption<bool>("nodes", "Build nodes datastructure"));
    add_option(new SimpleOption<bool>("edges", "Build edges datastructure"));
    add_option(new SimpleOption<bool>("cells", "Build cells datastructure"));
    add_option(new SimpleOption<bool>("brick", "Build brick dual mesh"));
    add_option(new SimpleOption<bool>("stats", "Write statistics file"));
    add_option(new SimpleOption<bool>("info", "Write Info"));
    add_option(new SimpleOption<bool>("periodic_x", "periodic mesh in x-direction"));
    add_option(new SimpleOption<bool>("periodic_y", "periodic mesh in y-direction"));
    add_option(new SimpleOption<bool>("torus", "Output mesh as torus"));
    add_option(new SimpleOption<bool>(
        "water", "Output elements containing water points (not specifying --water or --land enables both)"));
    add_option(new SimpleOption<bool>(
        "land", "Output elements containing land points (not specifying --water or --land enables both)"));
    add_option(new SimpleOption<bool>("fixup", "Apply custom fixes to the mesh where it applies"));
    add_option(new SimpleOption<bool>("gmsh", "Output gmsh (default=true)"));
    add_option(new SimpleOption<long>("partition", "partition [0:partitions]"));
    add_option(new SimpleOption<long>("partitions", "Number of partitions"));
    add_option(new SimpleOption<bool>("partition-graph", "Output partition graph"));
    add_option(new SimpleOption<bool>("partition-polygons", "Output partition polygons python visualization scripts"));
}

//-----------------------------------------------------------------------------

std::string get_arg(const AtlasTool::Args& args, const std::string& flag, const std::string& default_value = "") {
    for (int i = 0; i < args.count() - 1; ++i) {
        if (args(i) == flag) {
            return args(i + 1);
        }
    }
    if (not default_value.empty()) {
        return default_value;
    }
    throw_Exception("Could not find argument for flag " + flag);
}

int Meshgen2Gmsh::execute(const Args& args) {
    key = "";
    args.get("grid.name", key);

    nodes = false;
    args.get("nodes", nodes);
    edges = false;
    args.get("edges", edges);
    cells = false;
    args.get("cells", cells);
    stats = false;
    args.get("stats", stats);
    info = false;
    args.get("info", info);
    halo = 0;
    args.get("halo", halo);
    bool dim_3d = false;
    args.get("3d", dim_3d);
    brick = false;
    args.get("brick", brick);
    ghost = false;
    args.get("ghost", ghost);

    key      = get_positional_arg(args, 0);
    path_out = get_arg(args, "-o", args.count() > 1 ? get_positional_arg(args, 1) : "mesh.msh");


    //    if ( path_in_str.empty() && key.empty() ) {
    //        Log::warning() << "missing argument --grid.name or --grid.json" << std::endl;
    //        Log::warning() << "Usage: " << usage() << std::endl;
    //        return failed();
    //    }

    if (edges) {
        halo = std::max(halo, 1l);
    }

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
                Log::error() << "Could not generate mesh for grid defined in file \"" << path_in << "\"" << std::endl;
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

    Log::debug() << "Domain: " << grid.domain() << std::endl;
    if (auto g = StructuredGrid(grid)) {
        Log::debug() << "Periodic: " << g.periodic() << std::endl;
    }
    Log::debug() << "Spec: " << grid.spec() << std::endl;

    auto meshgenerator = make_meshgenerator(grid, args);
    auto partitioner   = make_partitioner(grid, args);

    Mesh mesh;
    try {
        Log::info() << "Generating mesh using " << meshgenerator.type() << " generator" << std::endl;
        mesh = meshgenerator.generate(grid, partitioner);
    }
    catch (eckit::Exception& e) {
        Log::error() << e.what() << std::endl;
        Log::error() << e.callStack() << std::endl;
        throw;
    }


    if ((grid.projection().units() == "degrees" && halo > 0) || nodes) {
        functionspace::NodeColumns nodes_fs(mesh, option::halo(halo));
    }
    else {
        Log::warning() << "Not yet implemented: building halo's with projections "
                          "not defined in degrees"
                       << std::endl;
        Log::warning() << "units: " << grid.projection().units() << std::endl;
    }
    if (brick) {
        build_brick_dual_mesh(grid, mesh);
    }
    if (edges && grid.projection().units() == "degrees") {
        functionspace::EdgeColumns edges_fs(mesh, option::halo(halo));
        if (not brick) {
            build_median_dual_mesh(mesh);
        }
    }

    if (cells) {
        functionspace::CellColumns cells_fs(mesh, option::halo(halo));
    }

    if (stats) {
        build_statistics(mesh);
    }

    bool lonlat             = args.getBool("lonlat", false);
    bool ij                 = args.getBool("ij", false);
    std::string coordinates = dim_3d ? "xyz" : lonlat ? "lonlat" : ij ? "ij" : "xy";
    args.get("coordinates", coordinates);

    if (args.getBool("gmsh", true)) {
        bool torus = false;
        args.get("torus", torus);
        if (torus) {
            if (auto g = StructuredGrid(grid)) {
                dim_3d = true;
                Log::debug() << "Building xyz representation for nodes on torus" << std::endl;
                mesh::actions::BuildTorusXYZField("xyz")(mesh, g.domain(), 5., 2., g.nxmax(), g.ny());
            }
            else {
                Log::error() << "Cannot output non StructuredGrid grids as torus at the moment" << std::endl;
            }
        }

        Config gmsh_config;
        gmsh_config.set("coordinates", coordinates);
        gmsh_config.set("edges", edges);
        gmsh_config.set("ghost", ghost);
        gmsh_config.set("info", info);
        if (args.has("land") || args.has("water")) {
            gmsh_config.set("land", args.getBool("land", false));
            gmsh_config.set("water", args.getBool("water", false));
        }
        atlas::output::Gmsh gmsh(path_out, gmsh_config);
        Log::info() << "Writing mesh to gmsh file \"" << path_out << "\" generated from grid \"" << grid.name() << "\""
                    << std::endl;
        gmsh.write(mesh);
    }

    if (info) {
        Log::info() << "Mesh partition footprint: " << eckit::Bytes(mesh.footprint()) << std::endl;
    }

    if (args.getBool("partition-graph", false)) {
        Log::info() << "Partitioning graph: \n" << mesh.partitionGraph() << std::endl;
    }

    if (args.getBool("partition-polygons", false)) {
        for (idx_t jhalo = 0; jhalo <= halo; ++jhalo) {
            mesh.polygon(jhalo).outputPythonScript("polygon_halo" + std::to_string(jhalo) + ".py",
                                                   Config("nodes", false)("coordinates", coordinates));
        }
    }

    return success();
}

//------------------------------------------------------------------------------

int main(int argc, char** argv) {
    Meshgen2Gmsh tool(argc, argv);
    return tool.start();
}
