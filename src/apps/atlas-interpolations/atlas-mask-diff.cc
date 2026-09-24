/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstddef>
#include <iomanip>
#include <string>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

#include "MaskIO.h"

namespace atlas {
namespace {

struct Labels {
    std::string first{"MASK1"};
    std::string second{"MASK2"};
    std::string zero{"off"};
    std::string one{"on"};
};

struct Comparison {
    std::size_t both_zero{0};
    std::size_t both_one{0};
    std::size_t only_first{0};
    std::size_t only_second{0};

    std::size_t equal() const { return both_zero + both_one; }
    std::size_t different() const { return only_first + only_second; }
};

enum class ComparisonValue : int {
    both_zero = 0,
    both_one = 1,
    only_first = 2,
    only_second = 3,
};

Labels read_labels(const AtlasTool::Args& args) {
    Labels labels;
    args.get("mask1-name", labels.first);
    args.get("mask2-name", labels.second);
    args.get("zero-name", labels.zero);
    args.get("one-name", labels.one);
    if (labels.first.empty() || labels.second.empty() || labels.zero.empty() || labels.one.empty()) {
        ATLAS_THROW_EXCEPTION("Mask and value names must not be empty");
    }
    return labels;
}

const std::string& value_name(int value, const Labels& labels) {
    return value == 0 ? labels.zero : labels.one;
}

Comparison compare(const std::vector<int>& first, const std::vector<int>& second) {
    Comparison result;
    for (std::size_t i = 0; i < first.size(); ++i) {
        if (first[i] == 0 && second[i] == 0) {
            ++result.both_zero;
        }
        else if (first[i] == 1 && second[i] == 1) {
            ++result.both_one;
        }
        else if (first[i] == 1) {
            ++result.only_first;
        }
        else {
            ++result.only_second;
        }
    }
    return result;
}

ComparisonValue comparison_value(int first, int second) {
    if (first == second) {
        return first == 0 ? ComparisonValue::both_zero : ComparisonValue::both_one;
    }
    return first == 1 ? ComparisonValue::only_first : ComparisonValue::only_second;
}

std::string comparison_field_name(const Labels& labels) {
    return "0=both_" + labels.zero + ", 1=both_" + labels.one + ", 2=" + labels.first + ":" + labels.one + "_" +
           labels.second + ":" + labels.zero + ", 3=" + labels.first + ":" + labels.zero + "_" + labels.second +
           ":" + labels.one;
}

void write_gmsh(const std::string& path, const MaskData& first, const MaskData& second, const Labels& labels) {
    Mesh mesh{Grid{first.grid}};
    functionspace::NodeColumns fs{mesh};
    Field field = fs.createField<int>(option::name(comparison_field_name(labels)));
    auto values = array::make_view<int, 1>(field);
    auto global_index = array::make_view<gidx_t, 1>(fs.global_index());
    auto ghost = array::make_view<int, 1>(fs.ghost());

    values.assign(-1);
    for (idx_t i = 0; i < values.size(); ++i) {
        const gidx_t index = global_index(i) - 1;
        if (index >= 0 && static_cast<std::size_t>(index) < first.values.size()) {
            values(i) = static_cast<int>(comparison_value(first.values[index], second.values[index]));
        }
        else if (!ghost(i)) {
            ATLAS_THROW_EXCEPTION("Mesh node has out-of-range global index " << global_index(i));
        }
    }
    fs.haloExchange(field);
    for (idx_t i = 0; i < values.size(); ++i) {
        if (values(i) < 0) {
            ATLAS_THROW_EXCEPTION("Could not map mesh node with global index " << global_index(i) << " to the grid");
        }
    }
    field.set_dirty(false);

    output::Gmsh gmsh{path, util::Config("coordinates", "xyz")};
    gmsh.write(mesh);
    gmsh.write(field);
}

}  // namespace

class AtlasMaskDiff : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Compare two Atlas mask records"; }
    std::string usage() override { return name() + " MASK1 MASK2 [OPTION]..."; }
    std::string longDescription() override {
        return "Exit status is 0 when the masks are equal and 1 when they differ.\n"
               "Gmsh values encode (MASK1,MASK2) as 0=(0,0), 1=(1,1), 2=(1,0), and 3=(0,1).";
    }
    int numberOfPositionalArguments() override { return 2; }
    int minimumPositionalArguments() override { return 2; }
    bool serial() override { return true; }

public:
    AtlasMaskDiff(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("mask1-name", "Display name for MASK1 (default: MASK1)"));
        add_option(new SimpleOption<std::string>("mask2-name", "Display name for MASK2 (default: MASK2)"));
        add_option(new SimpleOption<std::string>("zero-name", "Meaning of mask value 0 (default: off)"));
        add_option(new SimpleOption<std::string>("one-name", "Meaning of mask value 1 (default: on)"));
        add_option(new SimpleOption<bool>("show", "Show differing indices and coordinates"));
        add_option(new SimpleOption<long>("max-differences", "Maximum differences shown (default: 20; 0: unlimited)"));
        add_option(new SimpleOption<std::string>("output", "Write the comparison field to this Gmsh file"));
    }
};

int AtlasMaskDiff::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("atlas-mask-diff");

    const MaskData first = MaskIO::read(get_positional_arg(args, 0));
    const MaskData second = MaskIO::read(get_positional_arg(args, 1));
    const Labels labels = read_labels(args);

    if (first.grid != second.grid) {
        auto grid1 = atlas::Grid(first.grid);
        auto grid2 = atlas::Grid(second.grid);
        if (grid1.uid() != grid2.uid()) {
            ATLAS_THROW_EXCEPTION("Masks use different grids: '" << first.grid << "' and '" << second.grid << "'");
        }
    }
    if (first.values.size() != second.values.size()) {
        ATLAS_THROW_EXCEPTION("Masks have different sizes: " << first.values.size() << " and " << second.values.size());
    }

    const Comparison result = compare(first.values, second.values);
    Log::info() << "Grid: " << first.grid << "\n"
                << "Points: " << first.values.size() << "\n"
                << "Equal: " << result.equal() << " (both " << labels.zero << ": " << result.both_zero << ", both "
                << labels.one << ": " << result.both_one << ")\n"
                << "Different: " << result.different() << " (" << labels.first << "=" << labels.one << ", "
                << labels.second << "=" << labels.zero << ": " << result.only_first << "; " << labels.first << "="
                << labels.zero << ", " << labels.second << "=" << labels.one << ": " << result.only_second << ")"
                << std::endl;

    if (args.getBool("show", false) && result.different() != 0) {
        const Grid grid{first.grid};
        if (grid.size() != first.values.size()) {
            ATLAS_THROW_EXCEPTION("Grid '" << first.grid << "' has " << grid.size() << " points but masks contain "
                                            << first.values.size());
        }

        long maximum = 20;
        args.get("max-differences", maximum);
        if (maximum < 0) {
            ATLAS_THROW_EXCEPTION("--max-differences must be non-negative");
        }

        std::size_t shown = 0;
        auto point = grid.lonlat().begin();
        for (std::size_t i = 0; i < first.values.size(); ++i, ++point) {
            if (first.values[i] != second.values[i] && (maximum == 0 || shown < static_cast<std::size_t>(maximum))) {
                Log::info() << "index " << i << ", point " << std::setprecision(12) << *point << ": " << labels.first
                            << "=" << value_name(first.values[i], labels) << ", " << labels.second << "="
                            << value_name(second.values[i], labels) << std::endl;
                ++shown;
            }
        }
        if (shown < result.different()) {
            Log::info() << "... " << result.different() - shown
                        << " more differences (use --max-differences=0 to show all)" << std::endl;
        }
    }

    std::string output_path;
    if (args.get("output", output_path)) {
        Log::info() << "Writing comparison field to '" << output_path << "'" << std::endl;
        write_gmsh(output_path, first, second, labels);
    }
    return result.different() == 0 ? success() : failed();
}

}  // namespace atlas

int main(int argc, char** argv) {
    atlas::AtlasMaskDiff tool(argc, argv);
    return tool.start();
}
