/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "eckit/linalg/LinearAlgebra.h"

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/internals/AtlasTool.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/runtime/Log.h"

#include "PartitionedMesh.h"

using namespace atlas;


class AtlasParallelInterpolation : public AtlasTool {

    void execute(const AtlasTool::Args& args);
    std::string briefDescription() {
        return "Demonstration of parallel interpolation";
    }
    std::string usage() {
        return name() + " [--source-gridname=gridname] [--target-gridname=gridname] [OPTION]... [--help]";
    }

    int numberOfPositionalArguments() { return -1; }
    int minimumPositionalArguments() { return 0; }

public:

    AtlasParallelInterpolation(int argc, char* argv[]) : AtlasTool(argc, argv) {
        add_option(new SimpleOption<size_t>     ("log-rank",         "use specific MPI rank for logging (default 0)"));
        add_option(new SimpleOption<bool>       ("log-statistics",   "show simple statistics on source/target (default false)"));

        add_option(new SimpleOption<std::string>("method",           "interpolation method (default finite-element)"));
        add_option(new SimpleOption<std::string>("backend",          "linear algebra backend"));

        add_option(new SimpleOption<std::string>("source-gridname",  "source gridname"));
        add_option(new SimpleOption<std::string>("target-gridname",  "target gridname"));

        add_option(new SimpleOption<std::string>("source-mesh-partitioner",           "source mesh partitioner (equal_regions (default), ...)"));
        add_option(new SimpleOption<bool>       ("source-mesh-generator",             "source mesh generator (default structured)"));
        add_option(new SimpleOption<bool>       ("source-mesh-generator-triangulate", "source mesh generator triangulate option (default false)"));
        add_option(new SimpleOption<bool>       ("source-mesh-generator-angle",       "source mesh generator angle option (default false)"));
        add_option(new SimpleOption<size_t>     ("source-mesh-halo",                  "source mesh halo size (default 0)"));

        add_option(new SimpleOption<std::string>("target-mesh-partitioner",           "target mesh partitioner (PrePartitionedPolygon, PrePartitionedBruteForce)"));
        add_option(new SimpleOption<bool>       ("target-mesh-generator",             "target mesh generator (default structured)"));
        add_option(new SimpleOption<bool>       ("target-mesh-generator-triangulate", "target mesh generator triangulate option (default false)"));
        add_option(new SimpleOption<bool>       ("target-mesh-generator-angle",       "target mesh generator angle option (default false)"));

        add_option(new SimpleOption<size_t>     ("target-mesh-halo",                  "target mesh halo size (default 1)"));

        add_option(new SimpleOption<size_t>     ("k-nearest-neighbours",              "k nearest neighbours (default 1)"));


        add_option(new SimpleOption<bool>       ("polygons",                          "Output python script that plots partitions as polygons"));


        add_option(new SimpleOption<bool>       ("with-backward",                     "Also do backward interpolation (default false)"));

        add_option(new SimpleOption<bool>       ("fallback_to_2d",             "When the ray-tracing algorithm fails, we can either increase the halo of the source grid, or try to fallback to 2d."));

    }

};


void AtlasParallelInterpolation::execute(const AtlasTool::Args& args) {

    // Get user options
    std::string option;
    bool trigs = false;
    double angle = 0.;

    bool log_statistics = false;
    args.get("log-statistics", log_statistics);

    size_t log_rank = 0;
    args.get("log-rank", log_rank);

    if (eckit::mpi::comm().rank() != log_rank) {
        Log::reset();
    }

    if (args.get("backend", option)) {
        eckit::linalg::LinearAlgebra::backend(option);
    }

    size_t source_mesh_halo = 0;
    args.get("source-mesh-halo", source_mesh_halo);

    size_t target_mesh_halo = 1;
    args.get("target-mesh-halo", target_mesh_halo);
    Log::debug() << "target-mesh-halo " << target_mesh_halo << std::endl;


    // Generate and partition source & target mesh
    // source mesh is partitioned on its own, the target mesh uses (pre-partitioned) source mesh

    option = args.get("source-gridname", option)? option : "O16";
    atlas::grid::Grid src_grid(option);
    interpolation::PartitionedMesh src(
                args.get("source-mesh-partitioner",           option)? option : "equal_regions",
                args.get("source-mesh-generator",             option)? option : "structured",
                args.get("source-mesh-generator-triangulate", trigs)?  trigs  : false,
                args.get("source-mesh-generator-angle",       angle)?  angle  : 0. );


    option = args.get("target-gridname", option)? option : "O32";
    atlas::grid::Grid tgt_grid(option);
    interpolation::PartitionedMesh tgt(
                args.get("target-mesh-partitioner",           option)? option : "PrePartitionedPolygon",
                args.get("target-mesh-generator",             option)? option : "structured",
                args.get("target-mesh-generator-triangulate", trigs)?  trigs  : false,
                args.get("target-mesh-generator-angle",       angle)?  angle  : 0. );

    Log::info() << "Partitioning source grid" << std::endl;
    src.partition(src_grid);

    Log::info() << "Partitioning target grid" << std::endl;
    tgt.partition(tgt_grid, src);

    Log::info() << "Increasing source mesh halo by " << source_mesh_halo << " elements." << std::endl;
    functionspace::NodeColumns src_functionspace(src.mesh(), source_mesh_halo);
    Log::info() << "Increasing target mesh halo by " << target_mesh_halo << " elements." << std::endl;
    functionspace::NodeColumns tgt_functionspace(tgt.mesh(), target_mesh_halo);

    src.writeGmsh("src-mesh.msh");
    tgt.writeGmsh("tgt-mesh.msh");

    // Setup interpolator relating source & target meshes before setting a source FunctionSpace halo
    std::string interpolator_option = "finite-element";
    args.get("method", interpolator_option);


    Log::info() << "Computing forward interpolator" << std::endl;

    interpolation::method::Method::Ptr interpolator_forward(interpolation::method::MethodFactory::build(interpolator_option, args));

    interpolator_forward->setup(src.mesh(), tgt.mesh());

    bool with_backward = false;
    args.get("with-backward",with_backward);
    interpolation::method::Method::Ptr interpolator_backward;

    if( with_backward ) {
      Log::info() << "Computing backward interpolator" << std::endl;

      interpolator_backward.reset( interpolation::method::MethodFactory::build(interpolator_option, args));

      interpolator_backward->setup(tgt.mesh(), src.mesh());
    }


    // Create source FunctionSpace and fields

    field::FieldSet src_fields;
    {
        src_fields.add( src_functionspace.createField<double>("funny_scalar_1" /*, nb_levels=10*/) );
        src_fields.add( src_functionspace.createField<double>("funny_scalar_2" /*, nb_levels=10*/) );

        // Helper constants
        const double
                deg2rad = M_PI / 180.,
                c_lat = 0. * M_PI,
                c_lon = 1. * M_PI,
                c_rad = 2. * M_PI / 9.;

        atlas::array::ArrayView< double, 2 > lonlat( src.mesh().nodes().lonlat() );
        atlas::array::ArrayView< double, 1 >
                src_scalar_1(src_fields[0]),
                src_scalar_2(src_fields[1]);
        for (size_t j = 0; j < src.mesh().nodes().size(); ++j) {
            const double lon = deg2rad * lonlat(j, 0);  // (lon)
            const double lat = deg2rad * lonlat(j, 1);  // (lat)
            const double
                    c2 = std::cos(lat),
                    s1 = std::sin((lon-c_lon)/2.),
                    s2 = std::sin((lat-c_lat)/2.),
                    dist = 2.0 * std::sqrt( c2*s1*c2*s1 + s2*s2 );
            src_scalar_1(j) = dist < c_rad? 0.5 * (1. + std::cos(M_PI*dist/c_rad)) : 0.;
            src_scalar_2(j) = -src_scalar_1(j);


            double x = lonlat(j, 0) - 180.;
            double y = lonlat(j, 1);

            src_scalar_1(j) = -std::tanh(y/10.*std::cos(50./std::sqrt(x*x+y*y))-x/10.*std::sin(50./std::sqrt(x*x+y*y)));

        }
    }

    field::FieldSet tgt_fields;
    for (size_t i = 0; i < src_fields.size(); ++i) {
        tgt_fields.add( tgt_functionspace.createField<double>(src_fields[i].name()) );
    }

    src_functionspace.haloExchange(src_fields);

    Log::info() << "Writing input to interpolation to src.msh" << std::endl;
    src.writeGmsh("src.msh", &src_fields);


    Log::info() << "Interpolate forward" << std::endl;

    // interpolate (matrix-field vector multiply and  synchronize results) (FIXME: necessary?)
    interpolator_forward->execute(src_fields, tgt_fields);
    tgt_functionspace.haloExchange(tgt_fields);

    if( with_backward ) {
        Log::info() << "Interpolate backward" << std::endl;

        interpolator_backward->execute(tgt_fields, src_fields);
        src_functionspace.haloExchange(src_fields);
    }

    Log::info() << "Interpolations done" << std::endl;

    // Report simple statistics (on source & target)
    if (log_statistics) {
        double meanA, minA, maxA, meanB, minB, maxB;
        size_t nA, nB;

        for (size_t i = 0; i < src_fields.size(); ++i) {

            src_functionspace.minimum(src_fields[i], minA);
            src_functionspace.maximum(src_fields[i], maxA);
            src_functionspace.mean   (src_fields[i], meanA, nA);

            tgt_functionspace.minimum(tgt_fields[i], minB);
            tgt_functionspace.maximum(tgt_fields[i], maxB);
            tgt_functionspace.mean   (tgt_fields[i], meanB, nB);

            Log::info() << "Field '" << src_fields[i].name() << "' (N,min,mean,max):"
                        << "\n\tsource:\t" << nA << ",\t" << minA << ",\t" << meanA << ",\t" << maxA
                        << "\n\ttarget:\t" << nB << ",\t" << minB << ",\t" << meanB << ",\t" << maxB
                        << std::endl;
        }
    }


    // Output results
    Log::info() << "Writing forward interpolation results to tgt.msh" << std::endl;
    tgt.writeGmsh("tgt.msh", &tgt_fields);

    if( with_backward ) {
      Log::info() << "Writing backward interpolation results to src-back.msh" << std::endl;
      src.writeGmsh("src-back.msh", &src_fields);
    }
}


int main(int argc, char* argv[]) {
    AtlasParallelInterpolation tool(argc, argv);
    return tool.start();
}
