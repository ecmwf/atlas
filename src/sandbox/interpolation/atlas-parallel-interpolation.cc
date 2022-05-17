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

#include "eckit/config/Resource.h"
#include "eckit/log/Plural.h"

#include "PartitionedMesh.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/linalg/sparse/Backend.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/function/VortexRollup.h"

using namespace atlas;

class AtlasParallelInterpolation : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Demonstration of parallel interpolation"; }
    std::string usage() override {
        return name() +
               " [--source-gridname=gridname] "
               "[--target-gridname=gridname] [OPTION]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasParallelInterpolation(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<size_t>("log-rank", "use specific MPI rank for logging (default 0)"));
        add_option(new SimpleOption<bool>("log-statistics", "show simple statistics on source/target (default false)"));
        add_option(new SimpleOption<bool>("output-polygons", "Output Python script that plots partitions polygons"));

        add_option(new SimpleOption<std::string>("method", "interpolation method (default finite-element)"));
        add_option(
            new SimpleOption<std::string>("backward-method", "backward interpolation method (default finite-element)"));
        add_option(new SimpleOption<std::string>("backend", "linear algebra backend"));
        add_option(new SimpleOption<size_t>("k-nearest-neighbours", "k-nearest neighbours (default 1)"));
        add_option(new SimpleOption<bool>("with-backward", "Also do backward interpolation (default false)"));

        add_option(new SimpleOption<std::string>("source-gridname", "source gridname"));
        add_option(new SimpleOption<std::string>("source-mesh-partitioner",
                                                 "source mesh partitioner (equal_regions (default), ...)"));
        add_option(
            new SimpleOption<std::string>("source-mesh-generator", "source mesh generator (default structured)"));
        add_option(new SimpleOption<bool>("source-mesh-generator-triangulate",
                                          "source mesh generator triangulate option (default false)"));
        add_option(
            new SimpleOption<double>("source-mesh-generator-angle", "source mesh generator angle option (default 0.)"));
        add_option(new SimpleOption<size_t>("source-mesh-halo", "source mesh halo size (default 1)"));

        add_option(new SimpleOption<std::string>("target-gridname", "target gridname"));
        add_option(new SimpleOption<std::string>("target-mesh-partitioner",
                                                 "target mesh partitioner "
                                                 "(spherical-polygon, "
                                                 "lonlat-polygon, brute-force)"));
        add_option(new SimpleOption<bool>("target-mesh-generator", "target mesh generator (default structured)"));
        add_option(new SimpleOption<bool>("target-mesh-generator-triangulate",
                                          "target mesh generator triangulate option (default false)"));
        add_option(
            new SimpleOption<double>("target-mesh-generator-angle", "target mesh generator angle option (default 0.)"));
        add_option(new SimpleOption<size_t>("target-mesh-halo", "target mesh halo size (default 1)"));
        add_option(
            new SimpleOption<bool>("forward-interpolator-output", "Output forward interpolator's points and weights"));
        add_option(new SimpleOption<bool>("backward-interpolator-output",
                                          "Output backward interpolator's points and weights"));
        add_option(new SimpleOption<bool>("skip-halo-exchange", "Skip halo exchange"));
        add_option(new SimpleOption<double>("missing-value", "Missing value to be inserted when projection fails"));
    }
};

int AtlasParallelInterpolation::execute(const AtlasTool::Args& args) {
    // Get user options
    std::string option;
    bool trigs   = false;
    double angle = 0.;

    bool log_statistics = false;
    args.get("log-statistics", log_statistics);

    std::string interpolation_method = "finite-element";
    args.get("method", interpolation_method);


    if (args.get("backend", option)) {
        linalg::sparse::current_backend(option);
    }

    bool with_backward = false;
    args.get("with-backward", with_backward);

    // Generate and partition source & target mesh
    // source mesh is partitioned on its own, the target mesh uses
    // (pre-partitioned) source mesh

    auto source_gridname = args.getString("source-gridname", "O16");
    auto target_gridname = args.getString("target-gridname", "O32");
    Log::info() << "atlas-parallel-interpolation from source grid " << source_gridname << " to " << target_gridname
                << std::endl;
    Grid src_grid(source_gridname);

    idx_t source_mesh_halo = 0;
    args.get("source-mesh-halo", source_mesh_halo);

    interpolation::PartitionedMesh src(args.get("source-mesh-partitioner", option) ? option : "default",
                                       args.get("source-mesh-generator", option) ? option : "default",
                                       args.get("source-mesh-generator-triangulate", trigs) ? trigs : false,
                                       args.get("source-mesh-generator-angle", angle) ? angle : 0., true);

    Log::info() << "Partitioning source grid, halo of " << eckit::Plural(source_mesh_halo, "element") << std::endl;
    src.partition(src_grid);

    if (eckit::Resource<bool>("--output-polygons", false)) {
        src.mesh().polygon(0).outputPythonScript("src-polygons.py");
    }

    FunctionSpace src_functionspace;
    bool structured = false;
    for (auto& is_structured : {"structured-bicubic", "bicubic", "structured-bilinear", "bilinear"}) {
        if (interpolation_method == is_structured) {
            structured = true;
            break;
        }
    }
    if (structured) {
        src_functionspace =
            functionspace::StructuredColumns{src.mesh().grid(), option::halo(std::max<idx_t>(2, source_mesh_halo)) |
                                                                    util::Config("periodic_points", true)};
    }
    else {
        src_functionspace = functionspace::NodeColumns{src.mesh(), option::halo(source_mesh_halo)};
    }
    src.writeGmsh("src-mesh.msh");

    Grid tgt_grid(target_gridname);

    idx_t target_mesh_halo = args.getInt("target-mesh-halo", 0);

    interpolation::PartitionedMesh tgt(args.get("target-mesh-partitioner", option) ? option : "spherical-polygon",
                                       args.get("target-mesh-generator", option) ? option : "default",
                                       args.get("target-mesh-generator-triangulate", trigs) ? trigs : false,
                                       args.get("target-mesh-generator-angle", angle) ? angle : 0.,
                                       with_backward ? true : false);

    Log::info() << "Partitioning target grid, halo of " << eckit::Plural(target_mesh_halo, "element") << std::endl;
    tgt.partition(tgt_grid, src);

    functionspace::NodeColumns tgt_functionspace(tgt.mesh(), option::halo(target_mesh_halo));
    tgt.writeGmsh("tgt-mesh.msh");

    if (eckit::Resource<bool>("--output-polygons", false)) {
        tgt.mesh().polygon(0).outputPythonScript("tgt-polygons.py");
    }

    // Setup interpolator relating source & target meshes before setting a source
    // FunctionSpace halo
    Log::info() << "Computing forward/backward interpolator" << std::endl;

    Interpolation interpolator_forward(option::type(interpolation_method), src_functionspace, tgt_functionspace);
    Interpolation interpolator_backward;

    if (with_backward) {
        std::string backward_interpolation_method = "finite-element";
        args.get("method", backward_interpolation_method);
        Log::info() << "Computing backward interpolator" << std::endl;
        interpolator_backward =
            Interpolation(option::type(backward_interpolation_method), tgt_functionspace, src_functionspace);
    }

    if (args.getBool("forward-interpolator-output", false)) {
        interpolator_forward.print(Log::info());
    }

    // Create source FunctionSpace and fields

    FieldSet src_fields;
    {
        src_fields.add(src_functionspace.createField<double>(option::name("funny_scalar_1")));
        src_fields.add(src_functionspace.createField<double>(option::name("funny_scalar_2")));

        // Helper constants
        const double deg2rad = M_PI / 180., c_lat = 0. * M_PI, c_lon = 1. * M_PI, c_rad = 2. * M_PI / 9.;

        array::ArrayView<double, 2> lonlat = [&]() {
            if (auto fs = functionspace::NodeColumns(src_functionspace)) {
                return array::make_view<double, 2>(fs.nodes().lonlat());
            }
            else if (auto fs = functionspace::StructuredColumns(src_functionspace)) {
                return array::make_view<double, 2>(fs.xy());
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }();
        array::ArrayView<double, 1> src_scalar_1 = array::make_view<double, 1>(src_fields[0]),
                                    src_scalar_2 = array::make_view<double, 1>(src_fields[1]);

        ATLAS_ASSERT(src_scalar_1.shape(0) == lonlat.shape(0));
        for (idx_t j = 0; j < lonlat.shape(0); ++j) {
            const double lon = deg2rad * lonlat(j, 0);  // (lon)
            const double lat = deg2rad * lonlat(j, 1);  // (lat)
            const double c2 = std::cos(lat), s1 = std::sin((lon - c_lon) / 2.), s2 = std::sin((lat - c_lat) / 2.),
                         dist = 2.0 * std::sqrt(c2 * s1 * c2 * s1 + s2 * s2);
            src_scalar_1(j)   = dist < c_rad ? 0.5 * (1. + std::cos(M_PI * dist / c_rad)) : 0.;
            src_scalar_2(j)   = -src_scalar_1(j);

            //            double x = lonlat( j, 0 ) - 180.;
            //            double y = lonlat( j, 1 );

            //            src_scalar_1( j ) = -std::tanh( y / 10. * std::cos( 50. / std::sqrt( x * x + y * y ) ) -
            //                                            x / 10. * std::sin( 50. / std::sqrt( x * x + y * y ) ) );

            src_scalar_1(j) = util::function::vortex_rollup(lonlat(j, 0), lonlat(j, 1), 1.);
        }
    }

    FieldSet tgt_fields;
    for (idx_t i = 0; i < src_fields.size(); ++i) {
        auto tgt_field = tgt_fields.add(tgt_functionspace.createField<double>(option::name(src_fields[i].name())));
        double missing_value;
        if (args.get("missing-value", missing_value)) {
            tgt_field.metadata().set("missing_value", missing_value);
        }
    }

    if (args.getBool("skip-halo-exchange", false)) {
        src_fields.set_dirty(false);
    }
    else {
        src_functionspace.haloExchange(src_fields);
    }

    Log::info() << "Writing input to interpolation to src.msh" << std::endl;
    src.writeGmsh("src.msh", src_fields);

    Log::info() << "Interpolate forward" << std::endl;

    // interpolate (matrix-field vector multiply and  synchronize results) (FIXME:
    // necessary?)
    interpolator_forward.execute(src_fields, tgt_fields);

    if (args.getBool("skip-halo-exchange", false)) {
        tgt_fields.set_dirty(false);
    }
    else {
        tgt_functionspace.haloExchange(tgt_fields);
    }

    if (with_backward) {
        Log::info() << "Interpolate backward" << std::endl;

        interpolator_backward.execute(tgt_fields, src_fields);
        if (args.getBool("skip-halo-exchange", false)) {
            src_fields.set_dirty(false);
        }
        else {
            src_functionspace.haloExchange(src_fields);
        }
    }

    Log::info() << "Interpolations done" << std::endl;

    // Report simple statistics (on source & target)
    if (auto src_nodecolumns = functionspace::NodeColumns{src_functionspace}) {
        if (log_statistics) {
            double meanA, minA, maxA, meanB, minB, maxB;
            idx_t nA, nB;

            for (idx_t i = 0; i < src_fields.size(); ++i) {
                src_nodecolumns.minimum(src_fields[i], minA);
                src_nodecolumns.maximum(src_fields[i], maxA);
                src_nodecolumns.mean(src_fields[i], meanA, nA);
                tgt_functionspace.minimum(tgt_fields[i], minB);
                tgt_functionspace.maximum(tgt_fields[i], maxB);
                tgt_functionspace.mean(tgt_fields[i], meanB, nB);

                Log::info() << "Field '" << src_fields[i].name() << "' (N,min,mean,max):"
                            << "\n\tsource:\t" << nA << ",\t" << minA << ",\t" << meanA << ",\t" << maxA
                            << "\n\ttarget:\t" << nB << ",\t" << minB << ",\t" << meanB << ",\t" << maxB << std::endl;
            }
        }
    }

    // Output results
    Log::info() << "Writing forward interpolation results to tgt.msh" << std::endl;
    tgt.writeGmsh("tgt.msh", tgt_fields);

    if (with_backward) {
        Log::info() << "Writing backward interpolation results to src-back.msh" << std::endl;
        src.writeGmsh("src-back.msh", src_fields);
    }
    return success();
}

int main(int argc, char* argv[]) {
    AtlasParallelInterpolation tool(argc, argv);
    return tool.start();
}
