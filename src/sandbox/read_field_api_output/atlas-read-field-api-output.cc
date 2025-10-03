/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "HDF5Reader.h"
#include "transpositions.h"
#include "interpolate.h"
#include "output_gmsh.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"

namespace atlas {

class Program : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "brief description"; }
    std::string usage() override {
        return name() + " [OPTION] ... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    Program(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("file", "input file"));
        add_option(new SimpleOption<std::string>("dataset", "dataset"));
        add_option(new SimpleOption<long>("index", "dataset"));
        add_option(new SimpleOption<std::string>("rad.grid", "target grid to interpolate to"));
        add_option(new SimpleOption<std::string>("rad.nproma", "target grid to nproma interpolate to"));
        add_option(new SimpleOption<std::string>("gmsh.coordinates", "coordinates for gmsh output [lonlat, xyz]"));
        add_option(new SimpleOption<std::string>("interpolation", "coordinates for gmsh output [lonlat, xyz]"));
    }
};

util::Config get_subconfiguration(const eckit::Configuration& config, const std::string& key) {
    util::Config c;
    config.get(key, c);
    return c;
}

int Program::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");

    HDF5Reader hdf5_reader(args.getString("file"), args.getString("dataset"));
    util::Metadata metadata = hdf5_reader.read_metadata();
    Log::info() << util::Config{metadata}.json() << std::endl;

    std::string name = metadata.getString("name");
    int nproma       = metadata.getLong("nproma");
    int nflds        = metadata.getLong("nflds");
    DataType datatype(metadata.getString("datatype"));
    Grid IFS_grid(metadata.getString("grid"));

    functionspace::BlockStructuredColumns IFS_blocked_fs(IFS_grid, util::Config("nproma",nproma));
    
    Field IFS_blocked_f = IFS_blocked_fs.createField(option::name(name)|option::datatype(datatype)|option::levels(nflds));
    hdf5_reader.read_field(IFS_blocked_f);

    output_gmsh(IFS_blocked_f, "ifs.msh", get_subconfiguration(args, "gmsh"));

    // Now this needs to go to different grid with different nproma
    functionspace::BlockStructuredColumns rad_blocked_fs;
    Field rad_blocked_f;
    {
        auto rad_grid   = Grid(args.getString("rad.grid",IFS_grid.name()));
        auto rad_nproma = args.getLong("rad.nproma",nproma);
        if (rad_grid == IFS_grid && rad_nproma == nproma) {
            rad_blocked_fs = IFS_blocked_fs;
            rad_blocked_f = IFS_blocked_f;
        }
        else {
            rad_blocked_fs = functionspace::BlockStructuredColumns(rad_grid, grid::MatchingPartitioner(IFS_blocked_fs), util::Config("nproma",rad_nproma));
            rad_blocked_f = rad_blocked_fs.createField(IFS_blocked_f);
            if (rad_grid == IFS_grid) {
                copy_blocked_to_blocked(IFS_blocked_f, rad_blocked_f);
            }
            else {
                // Interpolation
                std::string interpolation_method = args.getString("interpolation","bicubic");
                interpolate(interpolation_method, IFS_blocked_f, rad_blocked_f);
            }
        }
    }

    output_gmsh(rad_blocked_f, "rad.msh", get_subconfiguration(args,"gmsh"));
    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::Program tool(argc, argv);
    return tool.start();
}
