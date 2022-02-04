/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"


using namespace atlas;


class AtlasParallelInterpolation : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Demonstration of parallel interpolation"; }
    std::string usage() override {
        return name() +
               " [--source=gridname] "
               "[--target=gridname] [OPTION]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasParallelInterpolation(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("source", "source gridname"));
        add_option(new SimpleOption<std::string>("target", "target gridname"));

        add_option(
            new SimpleOption<std::string>("method", "interpolation method [linear (default), cubic, quasicubic]"));
        add_option(new SimpleOption<bool>("with-matrix", "Store matrix for consecutive interpolations"));

        add_option(new SimpleOption<std::string>("partitioner", "source partitioner [equal_regions (default), ...]"));

        add_option(new SimpleOption<bool>("output-gmsh", "Output gmsh files src_field.msh and tgt_field.msh"));
        add_option(new SimpleOption<bool>(
            "output-polygons",
            "Output Python scripts src-polygons.py and tgt-polygons.py that plots partitions polygons"));
        add_option(
            new SimpleOption<bool>("output-inscribed-rectangle",
                                   "Output Python scripts src-inscribed-rectangle.py and tgt-inscribed-rectangle.py "
                                   "that plots iscribed rectangular domain for each partition"));

        add_option(
            new SimpleOption<std::string>("init", "Setup initial source field [ zero, vortex-rollup (default) ]"));
        add_option(new SimpleOption<long>("vortex-rollup", "Value that controls vortex rollup (default = 0)"));
        add_option(new SimpleOption<bool>("with-backwards", "Do backwards interpolation"));
    }
};

static Config processed_config(const eckit::Configuration& _config) {
    Config config;
    if (_config.has("partitioner")) {
        config.set("partitioner", option::type(_config.getString("partitioner")));
    }
    std::string scheme_str = _config.getString("method", "linear");
    if (scheme_str == "linear") {
        config.set("type", "structured-linear2D");
        config.set("halo", 1);
        // The stencil does not require any halo, but we set it to 1 for pole treatment!
    }
    if (scheme_str == "cubic") {
        config.set("type", "structured-cubic2D");
        config.set("halo", 2);
    }
    if (scheme_str == "quasicubic") {
        config.set("type", "structured-quasicubic2D");
        config.set("halo", 2);
    }
    config.set("name", scheme_str);
    config.set("matrix_free", not _config.getBool("with-matrix", false));
    return config;
}

int AtlasParallelInterpolation::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("AtlasParallelInterpolation::execute");
    auto source_gridname = args.getString("source", "O32");
    auto target_gridname = args.getString("target", "O64");
    auto require_polygon = args.getBool("output-polygons", false) || args.getBool("output-inscribed-rectangle", false);
    auto config          = processed_config(args);

    Log::info() << "atlas-parallel-interpolation from source grid " << source_gridname << " to " << target_gridname
                << std::endl;
    Grid src_grid(source_gridname);
    Grid tgt_grid(target_gridname);

    functionspace::StructuredColumns src_fs;
    functionspace::StructuredColumns tgt_fs;

    int nlev = 0;

    ATLAS_TRACE_SCOPE("Create source function space") {
        src_fs = functionspace::StructuredColumns{src_grid, config | option::levels(nlev)};
    }

    bool output_nodes = true;
    if (require_polygon) {
        const auto& src_poly = src_fs.polygon();
        if (args.getBool("output-polygons", false)) {
            src_poly.outputPythonScript("src-polygons.py", Config("nodes", output_nodes));
        }
        if (args.getBool("output-inscribed-rectangle", false)) {
            src_poly.outputPythonScript("src-inscribed-rectangle.py",
                                        Config("nodes", output_nodes) | Config("inscribed_rectangle", true));
        }
    }

    ATLAS_TRACE_SCOPE("Create target function space") {
        auto tgt_partitioner = grid::MatchingPartitioner(src_fs);
        tgt_fs = functionspace::StructuredColumns{tgt_grid, tgt_partitioner, config | option::levels(nlev)};
    }


    if (require_polygon) {
        const auto& tgt_poly = tgt_fs.polygon();
        if (args.getBool("output-polygons", false)) {
            tgt_poly.outputPythonScript("tgt-polygons.py", Config("nodes", output_nodes));
        }
        if (args.getBool("output-inscribed-rectangle", false)) {
            tgt_poly.outputPythonScript("tgt-inscribed-rectangle.py",
                                        Config("nodes", output_nodes)("inscribed_rectangle", true));
        }
    }


    Field src_field = src_fs.createField<double>(option::name("source"));
    Field tgt_field = tgt_fs.createField<double>(option::name("target"));

    // Initialize source
    {
        ATLAS_TRACE("Initialize source");
        auto lonlat = array::make_view<double, 2>(src_fs.xy());
        auto source = array::make_view<double, 1>(src_field);
        idx_t size  = src_fs.size();

        if (args.getString("init", "vortex-rollup") == "zero") {
            atlas_omp_parallel_for(idx_t n = 0; n < size; ++n) { source(n) = 0.; }
        }
        else {
            idx_t k = args.getInt("vortex-rollup", 0);
            atlas_omp_parallel_for(idx_t n = 0; n < size; ++n) {
                source(n) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
            }
        }
    }

    src_field.haloExchange();

    output::Output src_gmsh;
    output::Output tgt_gmsh;

    if (args.getBool("output-gmsh", false)) {
        src_gmsh = output::Gmsh("src_field.msh");
        src_gmsh.write(src_field);
    }

    {
        ATLAS_TRACE("Interpolation: Source to Target");
        Interpolation interpolation_fwd(config, src_fs, tgt_fs);
        interpolation_fwd.execute(src_field, tgt_field);
    }

    if (args.getBool("output-gmsh", false)) {
        tgt_field.haloExchange();
        tgt_gmsh = output::Gmsh("tgt_field.msh", Config("ghost", true));
        tgt_gmsh.write(tgt_field);
    }

    if (args.getBool("with-backwards", false)) {
        tgt_field.haloExchange();
        {
            ATLAS_TRACE("Interpolation: Target to Source");
            Interpolation interpolation_bwd(config, tgt_fs, src_fs);
            interpolation_bwd.execute(tgt_field, src_field);
        }

        if (args.getBool("output-gmsh", false)) {
            src_gmsh.write(src_field);
        }
    }

    ATLAS_TRACE_SCOPE("Load imbalance") { mpi::comm().barrier(); }

    return success();
}


int main(int argc, char* argv[]) {
    AtlasParallelInterpolation tool(argc, argv);
    return tool.start();
}
