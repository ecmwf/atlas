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
        add_option(new SimpleOption<std::string>("gmsh.coordinates", "coordinates for gmsh output [lonlat, xyz]"));
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
    int nflds        = metadata.getLong("nfld");
    DataType datatype(metadata.getString("datatype"));
    Grid IFS_grid(metadata.getString("grid"));

    functionspace::BlockStructuredColumns IFS_blocked_fs(IFS_grid, util::Config("nproma",nproma));
    
    Field IFS_blocked_f;
    if (nflds == 1) {
        IFS_blocked_f = IFS_blocked_fs.createField(option::name(name)|option::datatype(datatype));
    }
    else {
        IFS_blocked_f = IFS_blocked_fs.createField(option::name(name)|option::datatype(datatype)|option::levels(nflds));
    }
    hdf5_reader.read_field(IFS_blocked_f);
    output_gmsh(IFS_blocked_f, "ifs.msh", get_subconfiguration(args, "gmsh"));
    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::Program tool(argc, argv);
    return tool.start();
}
