/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <netcdf>

#include "atlas/grid.h"
#include "atlas/io/atlas-io.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/detail/filesystem.h"

#include "MaskIO.h"

namespace atlas {
namespace {

std::unique_ptr<netCDF::NcFile> open_netcdf(const std::string& path) {
    if (!filesystem::exists(path)) {
        ATLAS_THROW_EXCEPTION("Input file does not exist: '" << path << "'");
    }
    if (!filesystem::is_regular_file(path)) {
        ATLAS_THROW_EXCEPTION("Input path is not a regular file: '" << path << "'");
    }

    try {
        return std::make_unique<netCDF::NcFile>(path, netCDF::NcFile::read);
    }
    catch (const netCDF::exceptions::NcException& error) {
        ATLAS_THROW_EXCEPTION("Input file '" << path << "' is not a readable NetCDF file: " << error.what());
    }
}

netCDF::NcDim require_dimension(const netCDF::NcFile& file, const std::string& path, const std::string& name) {
    const netCDF::NcDim dimension = file.getDim(name);
    if (dimension.isNull()) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: missing dimension '" << name << "'");
    }
    if (dimension.getSize() == 0) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: dimension '" << name << "' is empty");
    }
    return dimension;
}

netCDF::NcVar require_variable(const netCDF::NcFile& file, const std::string& path, const std::string& name,
                              const std::vector<std::string>& dimensions) {
    const netCDF::NcVar variable = file.getVar(name);
    if (variable.isNull()) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: missing variable '" << name << "'");
    }
    if (variable.getDimCount() != static_cast<int>(dimensions.size())) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: variable '" << name << "' must have "
                                                << dimensions.size() << " dimension(s)");
    }
    for (std::size_t i = 0; i < dimensions.size(); ++i) {
        if (variable.getDim(static_cast<int>(i)).getName() != dimensions[i]) {
            ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: dimension " << i << " of variable '"
                                                    << name << "' must be '" << dimensions[i] << "'");
        }
    }
    return variable;
}

bool is_integer_type(nc_type type) {
    return type == NC_BYTE || type == NC_SHORT || type == NC_INT || type == NC_INT64 || type == NC_UBYTE ||
           type == NC_USHORT || type == NC_UINT || type == NC_UINT64;
}

struct ScripLayout {
    std::string side;
    std::string prefix;

    std::string name(const std::string& unprefixed) const { return prefix + unprefixed; }
};

bool has_layout(const netCDF::NcFile& file, const ScripLayout& layout) {
    return !file.getDim(layout.name("grid_size")).isNull() &&
           !file.getVar(layout.name("grid_imask")).isNull();
}

bool side_is_required(const netCDF::NcFile& file) {
    return file.getDim("grid_size").isNull() && has_layout(file, {"src", "src_"}) &&
           has_layout(file, {"dst", "dst_"});
}

ScripLayout select_layout(const netCDF::NcFile& file, const std::string& path, const std::string& side) {
    const ScripLayout standard{"", ""};
    const ScripLayout source{"src", "src_"};
    const ScripLayout destination{"dst", "dst_"};

    if (!side.empty()) {
        const ScripLayout layout = side == "src" ? source : side == "dst" ? destination : ScripLayout{};
        if (layout.side.empty()) {
            ATLAS_THROW_EXCEPTION("Invalid --side='" << side << "': expected 'src' or 'dst'");
        }
        if (!has_layout(file, layout)) {
            ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' does not contain the " << layout.side
                                                    << " SCRIP grid (expected '" << layout.name("grid_size")
                                                    << "' and '" << layout.name("grid_imask") << "')");
        }
        return layout;
    }

    if (has_layout(file, standard)) {
        return standard;
    }
    const bool has_source = has_layout(file, source);
    const bool has_destination = has_layout(file, destination);
    if (has_source && has_destination) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path
                                                << "' contains both source and destination SCRIP grids. "
                                                   "Select one with --side=src or --side=dst");
    }
    if (has_source) {
        return source;
    }
    if (has_destination) {
        return destination;
    }
    ATLAS_THROW_EXCEPTION("NetCDF file '" << path
                                            << "' contains no SCRIP grid. Expected 'grid_size', 'src_grid_size', "
                                               "or 'dst_grid_size' with a corresponding mask variable");
}

void validate_scrip(const netCDF::NcFile& file, const std::string& path, const ScripLayout& layout) {
    const std::string grid_size_name = layout.name("grid_size");
    const std::string grid_corners_name = layout.name("grid_corners");
    const std::string grid_rank_name = layout.name("grid_rank");
    const netCDF::NcDim grid_size = require_dimension(file, path, grid_size_name);
    require_dimension(file, path, grid_corners_name);
    const netCDF::NcDim grid_rank = require_dimension(file, path, grid_rank_name);

    const std::string grid_dims_name = layout.name("grid_dims");
    const netCDF::NcVar grid_dims = require_variable(file, path, grid_dims_name, {grid_rank_name});
    if (!is_integer_type(grid_dims.getType().getId())) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: variable '" << grid_dims_name
                                                << "' must have an integer type");
    }
    std::vector<long long> dimensions(grid_rank.getSize());
    grid_dims.getVar(dimensions.data());
    std::size_t represented_size = 1;
    for (long long dimension : dimensions) {
        if (dimension <= 0 || static_cast<std::size_t>(dimension) > grid_size.getSize() / represented_size) {
            ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: invalid value in '" << grid_dims_name
                                                    << "'");
        }
        represented_size *= static_cast<std::size_t>(dimension);
    }
    if (represented_size != grid_size.getSize()) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: product of '" << grid_dims_name << "' is "
                                                << represented_size << " but '" << grid_size_name << "' is "
                                                << grid_size.getSize());
    }

    for (const auto& variable :
         {require_variable(file, path, layout.name("grid_center_lat"), {grid_size_name}),
          require_variable(file, path, layout.name("grid_center_lon"), {grid_size_name}),
          require_variable(file, path, layout.name("grid_corner_lat"), {grid_size_name, grid_corners_name}),
          require_variable(file, path, layout.name("grid_corner_lon"), {grid_size_name, grid_corners_name})}) {
        const auto type = variable.getType().getId();
        if (type != NC_FLOAT && type != NC_DOUBLE) {
            ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: coordinate variable '"
                                                    << variable.getName() << "' must have a floating-point type");
        }
    }
    const std::string mask_variable = layout.name("grid_imask");
    const netCDF::NcVar mask = require_variable(file, path, mask_variable, {grid_size_name});
    if (!is_integer_type(mask.getType().getId())) {
        ATLAS_THROW_EXCEPTION("NetCDF file '" << path << "' is not SCRIP: mask variable '" << mask_variable
                                                << "' must have an integer type");
    }
}

std::string attribute_value(const netCDF::NcAtt& attribute) {
    std::ostringstream out;
    if (attribute.getType().getId() == NC_CHAR || attribute.getType().getId() == NC_STRING) {
        std::string value;
        attribute.getValues(value);
        out << '"' << value << '"';
    }
    else {
        std::vector<double> values(attribute.getAttLength());
        attribute.getValues(values.data());
        for (std::size_t i = 0; i < values.size(); ++i) {
            out << (i == 0 ? "" : ", ") << values[i];
        }
    }
    return out.str();
}

void list_netcdf_scrip(const std::string& path, const std::string& side) {
    const auto file = open_netcdf(path);

    Log::info() << "NetCDF file: " << path << "\nDimensions:" << std::endl;
    for (const auto& entry : file->getDims()) {
        const auto& dimension = entry.second;
        Log::info() << "  " << dimension.getName() << " = " << dimension.getSize()
                    << (dimension.isUnlimited() ? " (unlimited)" : "") << std::endl;
    }

    Log::info() << "Variables:" << std::endl;
    for (const auto& entry : file->getVars()) {
        const auto& variable = entry.second;
        Log::info() << "  " << variable.getType().getName() << " " << variable.getName() << "(";
        const auto dimensions = variable.getDims();
        for (std::size_t i = 0; i < dimensions.size(); ++i) {
            Log::info() << (i == 0 ? "" : ", ") << dimensions[i].getName() << "=" << dimensions[i].getSize();
        }
        Log::info() << ")" << std::endl;
        for (const auto& attribute : variable.getAtts()) {
            Log::info() << "    @" << attribute.first << " (" << attribute.second.getType().getName() << ") = "
                        << attribute_value(attribute.second) << std::endl;
        }
    }

    Log::info() << "Global attributes:" << std::endl;
    for (const auto& attribute : file->getAtts()) {
        Log::info() << "  @" << attribute.first << " (" << attribute.second.getType().getName() << ") = "
                    << attribute_value(attribute.second) << std::endl;
    }

    if (side.empty() && side_is_required(*file)) {
        for (const ScripLayout& layout : {ScripLayout{"src", "src_"}, ScripLayout{"dst", "dst_"}}) {
            validate_scrip(*file, path, layout);
            Log::info() << "SCRIP " << layout.side << " grid: valid\nMask variable: " << layout.name("grid_imask")
                        << std::endl;
        }
    }
    else {
        const ScripLayout layout = select_layout(*file, path, side);
        validate_scrip(*file, path, layout);
        Log::info() << "SCRIP conventions: valid\nMask variable: " << layout.name("grid_imask") << std::endl;
    }
}

MaskData extract_mask_scrip(const std::string& input_path, const std::string& grid_name, const std::string& side) {
    try {
        const auto file = open_netcdf(input_path);
        const ScripLayout layout = select_layout(*file, input_path, side);
        validate_scrip(*file, input_path, layout);
        const netCDF::NcDim grid_size = file->getDim(layout.name("grid_size"));

        const std::string mask_variable = layout.name("grid_imask");
        const netCDF::NcVar variable = file->getVar(mask_variable);

        MaskData data{grid_name, std::vector<int>(grid_size.getSize())};
        variable.getVar(data.values.data());
        MaskIO::validate(input_path + ":" + mask_variable, data.values);
        return data;
    }
    catch (const netCDF::exceptions::NcException& error) {
        ATLAS_THROW_EXCEPTION("Could not read SCRIP mask from '" << input_path << "': " << error.what());
    }
}

MaskData extract_mask_atlas_field(const std::string& input_path, long field_index, double missing_value) {
    if (field_index < 0) {
        ATLAS_THROW_EXCEPTION("--field-index must be non-negative");
    }
    if (!filesystem::exists(input_path)) {
        ATLAS_THROW_EXCEPTION("Input file does not exist: '" << input_path << "'");
    }

    std::string grid_name;
    std::vector<double> values;
    const std::string field = "fields[" + std::to_string(field_index) + "].array";
    try {
        io::RecordReader reader(input_path);
        reader.read("grid.name", grid_name);
        reader.read(field, values);
        reader.wait();
    }
    catch (const std::exception& error) {
        ATLAS_THROW_EXCEPTION("Could not read AtlasIO field '" << field << "' from '" << input_path
                                                                 << "': " << error.what());
    }

    const Grid grid{grid_name};
    if (values.size() != grid.size()) {
        ATLAS_THROW_EXCEPTION("AtlasIO field '" << field << "' contains " << values.size() << " values but grid '"
                                                  << grid_name << "' has " << grid.size() << " points");
    }

    MaskData data{grid_name, std::vector<int>(values.size())};
    std::transform(values.begin(), values.end(), data.values.begin(),
                   [missing_value](double value) { return value == missing_value ? 0 : 1; });
    return data;
}

}  // namespace

class AtlasMaskExtract : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Extract a mask into an Atlas mask record"; }
    std::string usage() override {
        return name() + " INPUT.nc OUTPUT.atlas GRID [--side=src|dst]\n" + name() +
               " INPUT.nc --list [--side=src|dst]\n" + name() +
             " INPUT.atlas OUTPUT.atlas --atlas-field [--field-index=N] [--missing-value=VALUE] (EXPERIMENTAL)";
    }
    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 1; }
    bool serial() override { return true; }

public:
    AtlasMaskExtract(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<bool>("list", "List NetCDF metadata and validate SCRIP conventions"));
        add_option(new SimpleOption<std::string>("side", "SCRIP remapping grid to extract: src or dst"));
        add_option(new SimpleOption<bool>("atlas-field", "EXPERIMENTAL: Extract a mask from an AtlasIO field"));
        add_option(new SimpleOption<long>("field-index", "AtlasIO field index (default: 0)"));
        add_option(new SimpleOption<double>("missing-value", "AtlasIO field missing value (default: 9999)"));
    }
};

int AtlasMaskExtract::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("atlas-mask-extract");

    const std::string input_path = get_positional_arg(args, 0);
    const std::string side = args.getString("side", "");
    const bool atlas_field = args.getBool("atlas-field", false);

    if (atlas_field) {
        Log::warning() << "\n"
                       << "**********************************************************************\n"
                       << "* WARNING: --atlas-field IS EXPERIMENTAL                             *\n"
                       << "* Its interface and mask extraction, semantics,                      *\n"
                       << "* and the file format itself may change.                             *\n"
                       << "* Verify the reported zero/one counts before using the output mask.  *\n"
                       << "**********************************************************************"
                       << std::endl;

        const std::string output_path = get_positional_arg(args, 1);
        if (output_path.empty() || !get_positional_arg(args, 2).empty()) {
            ATLAS_THROW_EXCEPTION("Usage: " << name()
                                              << " INPUT.atlas OUTPUT.atlas --atlas-field [--field-index=N] "
                                                 "[--missing-value=VALUE]");
        }
        if (args.getBool("list", false) || !side.empty()) {
            ATLAS_THROW_EXCEPTION("--atlas-field cannot be combined with --list or --side");
        }
        if (input_path == output_path) {
            ATLAS_THROW_EXCEPTION("Input and output paths must be different");
        }

        const long field_index = args.getLong("field-index", 0);
        const double missing_value = args.getDouble("missing-value", 9999.);
        const MaskData data = extract_mask_atlas_field(input_path, field_index, missing_value);
        MaskIO::write(output_path, data);

        const std::size_t ones = static_cast<std::size_t>(std::count(data.values.begin(), data.values.end(), 1));
        Log::info() << "Extracted " << data.values.size() << " values from fields[" << field_index
                    << "].array for grid " << data.grid << "\nZeros: " << data.values.size() - ones
                    << "\nOnes: " << ones << "\nWrote '" << output_path << "'" << std::endl;
        return success();
    }

    if (args.getBool("list", false)) {
        if (!get_positional_arg(args, 1).empty()) {
            ATLAS_THROW_EXCEPTION("Usage: " << name() << " INPUT.nc --list");
        }
        list_netcdf_scrip(input_path, side);
        return success();
    }

    if (side.empty() && side_is_required(*open_netcdf(input_path))) {
        Log::error() << "NetCDF file '" << input_path
                     << "' contains two valid SCRIP grids. Select the mask to extract with --side=src or --side=dst."
                     << std::endl;
        return failed();
    }

    const std::string output_path = get_positional_arg(args, 1);
    const std::string grid_name = get_positional_arg(args, 2);
    if (output_path.empty() || grid_name.empty() || !get_positional_arg(args, 3).empty()) {
        ATLAS_THROW_EXCEPTION("Usage: " << name() << " INPUT.nc OUTPUT.atlas GRID [--side=src|dst]");
    }
    if (input_path == output_path) {
        ATLAS_THROW_EXCEPTION("Input and output paths must be different");
    }

    const MaskData data = extract_mask_scrip(input_path, grid_name, side);
    const Grid grid{grid_name};
    if (grid.size() != data.values.size()) {
        ATLAS_THROW_EXCEPTION("Grid '" << grid_name << "' has " << grid.size() << " points but NetCDF mask contains "
                        << data.values.size() << " values");
    }
    MaskIO::write(output_path, data);

    const std::size_t ones = static_cast<std::size_t>(std::count(data.values.begin(), data.values.end(), 1));
    Log::info() << "Extracted " << data.values.size() << " values for grid " << grid_name << "\nZeros: "
                << data.values.size() - ones << "\nOnes: " << ones << "\nWrote '" << output_path << "'" << std::endl;
    return success();
}

}  // namespace atlas

int main(int argc, char** argv) {
    atlas::AtlasMaskExtract tool(argc, argv);
    return tool.start();
}
