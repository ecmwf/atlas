/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "HDF5Reader.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"

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
    }

};

int Program::execute(const AtlasTool::Args& args) {
    std::string input_file   = args.getString("file");
    std::string dataset_name = args.getString("dataset");
    int field_index = int(args.getLong("index", 0));

    // Open the HDF5 file
    HDF5_File h5_file(input_file);
    if (!h5_file.is_open()) {
        return failed();
    }

    // Get the dataset from the file
    HDF5_Dataset h5_dataset = h5_file.dataset(dataset_name);
    if (!h5_dataset.is_open()) {
        return failed();
    }
    Log::info() << "Dataset name     : " << h5_dataset.name() << '\n';
    Log::info() << "Dataset rank     : " << h5_dataset.rank() << '\n';
    Log::info() << "Dataset extents  : " << h5_dataset.extents() << '\n';
    Log::info() << "Dataset datatype : " << h5_dataset.datatype().str() << '\n';

    // Create and read in field
    Field field(h5_dataset.name(), h5_dataset.datatype(), h5_dataset.extents());
    h5_dataset.read(field);

    Log::info() << "field.name     : " << field.name() << '\n';
    Log::info() << "field.rank     : " << field.rank() << '\n';
    Log::info() << "field.shape    : " << field.shape() << '\n';
    Log::info() << "field.datatype : " << field.datatype().str() << '\n';

    int nproma = field.shape(field.rank()-1);
    int nblk   = field.shape(0);
    int nfld   = field.rank() == 3 ? field.shape(1) : 1;

    Log::info() << "nproma : " << nproma << '\n';
    Log::info() << "nblk   : " << nblk << '\n';
    Log::info() << "nfld   : " << nfld << '\n';
    Log::info() << "nblk * nproma : " << nblk * nproma << '\n';

    // Print
    if (false)
    {
        if (field.rank() == 1 && field.datatype().kind() == DataType::KIND_REAL64) {
            array::make_view<double,1>(field).dump(Log::info());
        }
        if (field.rank() == 2 && field.datatype().kind() == DataType::KIND_REAL64) {
            array::make_view<double,2>(field).dump(Log::info());
        }
        if (field.rank() == 3 && field.datatype().kind() == DataType::KIND_REAL64) {
            array::make_view<double,3>(field).slice(0,array::Range::all(),array::Range::all()).dump(Log::info());
        }
    }

    std::string gridname;
    if (nproma * nblk == 2048) {
        gridname = "F16";
    }
    if (nproma * nblk == 35720) {
        gridname = "N80"; // Grid has actually size 35718, but nproma*nblk contains padding
    }
    if (gridname.empty()) {
        Log::error() << "ERROR: Could not detect grid" << std::endl;
        return failed();
    }
    Log::info() << "gridname = " << gridname << std::endl;
    auto grid = Grid(gridname);
    auto mesh = Mesh(grid);
    auto fs = functionspace::BlockStructuredColumns(grid, util::Config("nproma",nproma));
    ATLAS_ASSERT(fs.nproma() == nproma);
    ATLAS_ASSERT(fs.nblks() == nblk);

    Log::info() << "Field for a single index:" << std::endl;
    auto f   = fs.createField(option::name("f")  | option::datatype(field.datatype()));
    auto fg  = fs.createField(option::name("fg") | option::datatype(field.datatype()) | option::global());
    Log::info() << "  f.shape() = " << f.shape() << std::endl;
    // Log::info() << "  fg.shape() = " << fg.shape() << std::endl;

    auto make_view = [](Field field, int field_index) {
        if (field.rank() == 3) {
            return array::make_view<double,3>(field).slice(array::Range::all(), field_index, array::Range::all());
        }
        else {
            return array::make_view<double,2>(field).slice(array::Range::all(),array::Range::all());
        }
    };

    if (true) {
        auto fieldv = make_view(field, field_index);
        auto fv     = array::make_view<double,2>(f);
        fv.assign(fieldv);
        fs.gather(f, fg);
    }
    else {
        auto fieldv = make_view(field, field_index);
        auto fgv    = array::make_view<double,1>(fg);
        int nmax = fg.size();
        for (int jblk=0; jblk<nblk; ++jblk) {
            for (int jrof=0; jrof<nproma; ++jrof) {
                int n = jblk*nproma + jrof;
                if (n < nmax) {
                    fgv(n) = fieldv(jblk, jrof);
                }
            }
        }
    }

    auto gmsh = output::Gmsh("mesh.msh");
    gmsh.write(mesh);
    gmsh.write(fg, functionspace::StructuredColumns(grid));
    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::Program tool(argc, argv);
    return tool.start();
}
