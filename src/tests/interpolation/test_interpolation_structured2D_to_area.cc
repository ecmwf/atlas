/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/PolygonXY.h"

#include "atlas/domain/Domain.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid/Partitioner.h"

#include "atlas/util/GridPointsJSONWriter.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::string input_gridname(const std::string& default_grid) {
    return eckit::Resource<std::string>("--input-grid", default_grid);
}

int output_rank(const int& default_rank) {
    return eckit::Resource<int>("--output-rank", default_rank);
}

std::string output_functionspace(const std::string& default_fs) {
    return eckit::Resource<std::string>("--output-functionspace", default_fs);
}


FieldSet create_source_fields(functionspace::StructuredColumns& fs, idx_t nb_fields, idx_t nb_levels) {
    using Value = double;
    FieldSet fields_source;
    auto lonlat = array::make_view<double, 2>(fs.xy());
    for (idx_t f = 0; f < nb_fields; ++f) {
        auto field_source = fields_source.add(fs.createField<Value>());
        auto source       = array::make_view<Value, 2>(field_source);
        for (idx_t n = 0; n < fs.size(); ++n) {
            for (idx_t k = 0; k < nb_levels; ++k) {
                source(n, k) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
            }
        };
    }
    return fields_source;
}
FieldSet create_target_fields(FunctionSpace& fs, idx_t nb_fields, idx_t nb_levels) {
    using Value = double;
    FieldSet fields_target;
    for (idx_t f = 0; f < nb_fields; ++f) {
        fields_target.add(fs.createField<Value>(option::levels(nb_levels)));
    }
    return fields_target;
}

void do_test( std::string type, int input_halo, bool expect_failure ) {
    idx_t outrank = output_rank(-1);
    idx_t nb_fields = 2;
    idx_t nb_levels = 19;

    // Setup Grid and functionspace
    Grid inputGrid(input_gridname("O40"));
    functionspace::StructuredColumns inputFS(inputGrid, option::levels(nb_levels) | option::halo(input_halo));

    inputFS.polygon(0).outputPythonScript("input.py");

    // Setup source field_set
    FieldSet fields_source = create_source_fields(inputFS, nb_fields, nb_levels);

    // Define cutout area and grid
    double boundingBoxNorth = 45;
    double boundingBoxSouth = -45;
    double boundingBoxEast = 140;
    double boundingBoxWest = 50;

    RectangularLonLatDomain rd({boundingBoxWest, boundingBoxEast}, {boundingBoxSouth, boundingBoxNorth});
    Grid areaGrid(inputGrid, rd);

    util::GridPointsJSONWriter input_writer(inputGrid,
        util::Config
            ("partition",mpi::rank())
            ("partitioner.partitions",mpi::size())
            ("field","lonlat")
    );
    if( mpi::rank() == outrank ) {
        Log::info() << "input grid coordinates of rank " << outrank << ": \n";
        input_writer.write(Log::info());
    }

    util::GridPointsJSONWriter output_writer(areaGrid, grid::MatchingPartitioner(inputFS),
        util::Config
            ("partition",mpi::rank())
            ("field","lonlat")
    );
    if( mpi::rank() == outrank ) {
        Log::info() << "output grid coordinates of rank " << outrank << ": \n";
        output_writer.write(Log::info());
    }

    FunctionSpace outputFS;
    std::string output_functionspace_type = output_functionspace("NodeColumns");

    if (output_functionspace_type=="NodeColumns") {
        Mesh areaMesh = MeshGenerator("structured").generate( areaGrid, grid::MatchingPartitioner(inputFS) );
        output::Gmsh gmsh{"area.msh"};
        gmsh.write(areaMesh);
        outputFS = atlas::functionspace::NodeColumns{areaMesh};
    }
    else if (output_functionspace_type=="StructuredColumns") {
        outputFS = functionspace::StructuredColumns(areaGrid, grid::MatchingPartitioner(inputFS) );
    }
    else if (output_functionspace_type=="PointCloud") {
        outputFS = atlas::functionspace::PointCloud{areaGrid, grid::MatchingPartitioner(inputFS)};
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    // setup interpolation
    Interpolation interpolation(option::type(type)|Config("matrix_free",false), inputFS, outputFS);

    // setup target field_set
    FieldSet fields_target = create_target_fields(outputFS, nb_fields, nb_levels);

    // execute interpolation
    interpolation.execute(fields_source, fields_target);

    if (atlas::functionspace::NodeColumns(outputFS)) {
        output::Gmsh gmsh{"area.msh","a"};
        gmsh.write(fields_target);
    }
}

CASE("test structured-linear2D, halo 3") {
    EXPECT_NO_THROW( do_test("structured-linear2D", 3, false) );
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
