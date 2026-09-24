/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

#include "MaskIO.h"

namespace atlas {
namespace {

struct InterpolationOptions {
    std::string type;
    long neighbours{0};

    util::Config config() const {
        util::Config result("type", type);
        if (type == "k-nearest-neighbours") {
            result.set(type, neighbours);
        }
        return result;
    }
};

std::string atlas_output_path(std::string path) {
    const std::string extension{".atlas"};
    if (path.size() < extension.size() || path.substr(path.size() - extension.size()) != extension) {
        ATLAS_THROW_EXCEPTION("Output mask must use the '.atlas' extension: " << path);
    }
    return path;
}

std::string resolve_output_path(const AtlasTool::Args& args, const Grid& target_grid) {
    std::string long_option;
    const bool has_long_option = args.get("output", long_option);

    std::string short_option;
    for (int i = 0; i < args.count(); ++i) {
        if (!args(i).empty() && args(i)[0] == '-' && args(i) != "-o") {
            ATLAS_THROW_EXCEPTION("Unknown option '" << args(i) << "'");
        }
        if (args(i) == "-o") {
            if (!short_option.empty()) {
                ATLAS_THROW_EXCEPTION("-o may only be specified once");
            }
            if (i + 1 >= args.count() || args(i + 1).empty() || args(i + 1)[0] == '-') {
                ATLAS_THROW_EXCEPTION("-o requires an output path");
            }
            short_option = args(i + 1);
        }
    }

    if (has_long_option && !short_option.empty()) {
        ATLAS_THROW_EXCEPTION("Use either -o or --output, not both");
    }
    if (has_long_option) {
        return atlas_output_path(long_option);
    }
    if (!short_option.empty()) {
        return atlas_output_path(short_option);
    }
    return "mask_" + target_grid.name() + ".atlas";
}

InterpolationOptions interpolation_options(const AtlasTool::Args& args) {
    InterpolationOptions options{"nearest-neighbour"};
    args.get("type", options.type);
    if (options.type.empty()) {
        ATLAS_THROW_EXCEPTION("Interpolation type must not be empty");
    }

    if (options.type == "nn" || options.type == "nearest") {
        options.type = "nearest-neighbour";
    }
    else if (options.type == "knn" || options.type == "k-nearest") {
        options.type = "k-nearest-neighbours";
    }

    if (options.type == "k-nearest-neighbours") {
        options.neighbours = 1;
        args.get("knn", options.neighbours);
        if (options.neighbours < 1) {
            ATLAS_THROW_EXCEPTION("--knn must be positive");
        }
    }

    return options;
}

}  // namespace

class AtlasMaskInterpolate : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Interpolate an Atlas mask to another grid"; }
    std::string usage() override { return name() + " SOURCE.atlas TARGET_GRID [-o OUTPUT.atlas] [OPTION]..."; }
    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 2; }

public:
    AtlasMaskInterpolate(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>(
            "type", "Interpolation type: nearest-neighbour (nn, nearest) or k-nearest-neighbours (knn, k-nearest)"));
        add_option(new SimpleOption<long>("knn", "Number of neighbours for --type=knn (default: 1)"));
        add_option(new SimpleOption<std::string>("output", "Output mask path (default: mask_<TARGET_GRID>.atlas)"));
        add_option(new SimpleOption<bool>("output-gmsh", "Write the interpolated field as Gmsh"));
    }
};

int AtlasMaskInterpolate::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("atlas-mask-interpolate");

    const std::string source_path = get_positional_arg(args, 0);
    const MaskData source_mask = MaskIO::read(source_path);
    const Grid source_grid{source_mask.grid};
    const Grid target_grid{get_positional_arg(args, 1)};
    if (!get_positional_arg(args, 2).empty()) {
        ATLAS_THROW_EXCEPTION("Usage: " << usage());
    }
    const std::string output_path = resolve_output_path(args, target_grid);
    const bool output_gmsh = args.getBool("output-gmsh", false);
    const std::string gmsh_path =
        output_gmsh ? output_path.substr(0, output_path.size() - std::string(".atlas").size()) + ".msh" : "";
    const InterpolationOptions interpolation = interpolation_options(args);
    if (source_grid.size() != source_mask.values.size()) {
        ATLAS_THROW_EXCEPTION("Grid '" << source_mask.grid << "' has " << source_grid.size()
                                        << " points but source mask contains " << source_mask.values.size());
    }

    if (mpi::rank() == 0) {
        Log::info() << "Interpolating mask\n"
                    << "  Source: " << source_path << "\n"
                    << "  Source grid: " << source_grid.name() << " (" << source_grid.size() << " points)\n"
                    << "  Target grid: " << target_grid.name() << " (" << target_grid.size() << " points)\n"
                    << "  Method: " << interpolation.type;
        if (interpolation.type == "k-nearest-neighbours") {
            Log::info() << " (neighbours: " << interpolation.neighbours << ")";
        }
        Log::info() << "\n  Output: " << output_path;
        if (output_gmsh) {
            Log::info() << "\n  Gmsh output: " << gmsh_path;
        }
        Log::info() << std::endl;
    }

    Interpolation interpolator(interpolation.config(), source_grid, target_grid);
    auto source_field = interpolator.source().createField<double>(option::name("mask"));
    auto source_values = array::make_view<double, 1>(source_field);
    auto source_global_index = array::make_view<gidx_t, 1>(interpolator.source().global_index());
    for (idx_t i = 0; i < source_values.size(); ++i) {
        const gidx_t index = source_global_index(i) - 1;
        if (index < 0 || static_cast<std::size_t>(index) >= source_mask.values.size()) {
            ATLAS_THROW_EXCEPTION("Source function-space point has out-of-range global index " << source_global_index(i));
        }
        source_values(i) = source_mask.values[index];
    }

    auto target_field = interpolator.target().createField<double>(option::name("mask"));
    interpolator.execute(source_field, target_field);

    if (output_gmsh) {
        Log::info() << "Writing interpolated field to '" << gmsh_path << "'" << std::endl;
        output::Gmsh gmsh(gmsh_path, util::Config("coordinates", "lonlat") | util::Config("ghost", false));
        if (functionspace::NodeColumns nodes = interpolator.target()) {
            gmsh.write(nodes.mesh());
        }
        gmsh.write(target_field);
    }

    auto global_field = interpolator.target().createField(target_field, option::global());
    interpolator.target().gather(target_field, global_field);
    if (mpi::rank() == 0) {
        auto global_values = array::make_view<double, 1>(global_field);
        MaskData target_mask{target_grid.name(), std::vector<int>(global_values.size())};
        for (idx_t i = 0; i < global_values.size(); ++i) {
            target_mask.values[i] = static_cast<int>(std::lround(global_values(i)));
        }
        MaskIO::write(output_path, target_mask);
        Log::info() << "Wrote " << target_mask.values.size() << " values for grid " << target_mask.grid << " to '"
                    << output_path << "'" << std::endl;
    }
    mpi::comm().barrier();
    return success();
}

}  // namespace atlas

int main(int argc, char** argv) {
    atlas::AtlasMaskInterpolate tool(argc, argv);
    return tool.start();
}
