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
        add_option(new SimpleOption<std::string>("variable", "dataset"));
    }

};

int Program::execute(const AtlasTool::Args& args) {
    std::string input_file   = args.getString("file");
    std::string dataset_name = args.getString("variable");

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

    // Print
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

    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::Program tool(argc, argv);
    return tool.start();
}
