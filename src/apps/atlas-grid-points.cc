/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>

#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/GridPointsJSONWriter.h"

//------------------------------------------------------------------------------

using namespace atlas;
using eckit::PathName;

//------------------------------------------------------------------------------

class Program : public AtlasTool {
    int execute(const Args& args) override;
    std::string briefDescription() override { return "Write grid points to file"; }
    std::string usage() override { return name() + " GRID [OPTION]... [--help]"; }
    std::string longDescription() override {
        return "    The 'GRID' argument can be either the name of a named grid, or the path to a"
               " YAML configuration file that describes the grid.\n"
               "Example values for grid names are: N80, F40, O24, L64x33, CS-ED-12. See the program "
               "'atlas-grids' for a list of named grids.\n"
               "\n";
    }

public:
    Program(int argc, char** argv);

private:
};

//-----------------------------------------------------------------------------

Program::Program(int argc, char** argv): AtlasTool(argc, argv) {
    add_option(new SimpleOption<std::string>("output.file", "Output file. If not specified, output is directed to stdout"));
    add_option(new SimpleOption<std::string>("output.format", "Output format. If not specified: json"));
    add_option(new SimpleOption<long>("field_base", "Base used for field output. Default=0"));
    add_option(new SimpleOption<long>("index_base", "Base used for index input. Default=0"));
    add_option(new SimpleOption<std::string>("index", 
        "Select grid point indices (first_index=<index_base>, last_index = <index_base> + <size> - 1). "
        "If not provided, all points are selected. "
        "Format: comma separated list, where '-' can be used to represent a range. e.g. '[1-3,5,7-10]'."
        "Square brackets are optional. white-spaces and newline characters are allowed as in a valid JSON array."));
    add_option(new SimpleOption<std::string>("field","Field to output. [\"lonlat\"(D),\"index\",\"partition\"]"));
    add_option(new Separator("Advanced"));
    add_option(new SimpleOption<std::string>("partitioner.type",
        "Partitioner [equal_regions,equal_area,checkerboard,equal_bands,regular_bands]"));
    add_option(new SimpleOption<long>("partition", "partition [0:partitions-1]"));
    add_option(new SimpleOption<long>("partitions", "Number of partitions"));
    add_option(new SimpleOption<long>("json.precision", "Number of digits after decimal in output"));
    add_option(new SimpleOption<bool>("json.pretty", "Pretty printing of json output"));
    add_option(new SimpleOption<bool>("verbose", "Output progress to stdout, default=false"));
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

int Program::execute(const Args& args) {

    Grid grid;
    {
        std::string key = args.count() ? args(0) : "";
        if (!key.empty()) {
            eckit::PathName path{key};
            grid = path.exists() ? Grid(Grid::Spec{path}) : Grid(key);

        }
    }
    if (!grid) {
        Log::error() << "Grid not specified as positional argument" << std::endl;
        return failed();
    }

    util::GridPointsJSONWriter writer{grid,args};
    if (mpi::rank() == 0 ) {
        std::string output_file;
        if (args.get("output.file",output_file)) {
            Log::info() << "Grid contains " << grid.size() << " points." << std::endl;
            std::ofstream out(output_file);
            writer.write(out, Log::info());
            Log::info() << "File " << output_file << " written." << std::endl;
        }
        else {
            writer.write(std::cout);
        }
    }
    return success();
}

//------------------------------------------------------------------------------

int main(int argc, char** argv) {
    Program tool(argc, argv);
    return tool.start();
}
