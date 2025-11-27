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
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/PolygonXY.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

static Grid get_input_grid() {
    Config config;
    config.set("type", "regional");
    config.set("lonlat(xmin,ymin)", std::vector<double>{4.03896, 51.1994});
    config.set("nx", 49);
    config.set("ny", 69);
    config.set("dx", 2500.);
    config.set("dy", 2500.);
    config.set("projection.type", "lambert_conformal_conic");
    config.set("projection.latitude0", 51.967);
    config.set("projection.longitude0", 4.9);
    return Grid{config};
/*

This grid extracted from a GRIB2 file received by Ulf Andrae to investigate support of LAM grids in the IFS using Atlas

***** FILE: grib_example_lambert.grib2
#==============   MESSAGE 1 ( length=208 )                 ==============
GRIB {
  # Meteorological products (grib2/tables/15/0.0.table)
  discipline = 0;
  editionNumber = 2;
  # French Weather Service - Toulouse (common/c-11.table)
  centre = 85;
  subCentre = 84;
  # Start of forecast (grib2/tables/15/1.2.table)
  significanceOfReferenceTime = 1;
  dataDate = 20251104;
  dataTime = 0;
  # Operational products (grib2/tables/15/1.3.table)
  productionStatusOfProcessedData = 0;
  # Forecast products (grib2/tables/15/1.4.table)
  typeOfProcessedData = 1;
  numberOfDataPoints = 3381;
  # There is no appended list (grib2/tables/15/3.11.table)
  interpretationOfNumberOfPoints = 0;
  # Unknown code table entry (grib2/tables/15/3.1.table)
  gridDefinitionTemplateNumber = 33;
  # Earth assumed spherical with radius of 6 371 229.0 m (grib2/tables/15/3.2.table)
  shapeOfTheEarth = 6;
  Nx = 49;
  Ny = 69;
  latitudeOfFirstGridPointInDegrees = 51.1994;
  longitudeOfFirstGridPointInDegrees = 4.03896;
  LaDInDegrees = 51.967;
  LoVInDegrees = 4.9;
  DxInMetres = 2500;
  DyInMetres = 2500;
  # (1=0)  North Pole is on the projection plane ;(2=0)  Only one projection centre is used :grib2/tables/15/3.5.table
  # flags: 00000000
  projectionCentreFlag = 0;
  iScansNegatively = 0;
  jScansPositively = 1;
  jPointsAreConsecutive = 0;
  alternativeRowScanning = 0;
  Latin1InDegrees = 51.967;
  Latin2 = 51967000;
  Latin2InDegrees = 51.967;
  latitudeOfSouthernPoleInDegrees = -90;
  longitudeOfSouthernPoleInDegrees = 0;
  Nux = 49;
  Ncx = 8;
  Nuy = 69;
  Ncy = 8;
  gridType = lambert_lam;
  NV = 0;
*/
}

static Grid get_output_grid(const std::string& name) {
    if (name == "lonlat_centred") {
        Config config;
        config.set("type", "regional");
        config.set("lonlat(xmin,ymin)", std::vector<double>{4.9-0.5, 51.967-0.5});
        config.set("nx", 41);
        config.set("ny", 41);
        config.set("dx", 0.5/20.);
        config.set("dy", 0.5/20.);
        config.set("projection.type", "lonlat");
        return Grid{config};
    }
    else if (name == "lonlat_bottomleft") {
        Config config;
        config.set("type", "regional");
        config.set("lonlat(xmin,ymin)", std::vector<double>{4.2, 51.3});
        config.set("nx", 21);
        config.set("ny", 21);
        config.set("dx", 0.5/20.);
        config.set("dy", 0.5/20.);
        config.set("projection.type", "lonlat");
        return Grid{config};
    }
    ATLAS_NOTIMPLEMENTED;
}

void check_input_structured_partition_polygon(FunctionSpace input_fs) {
    ATLAS_TRACE();
    input_fs.polygon().outputPythonScript("input_partitions_n"+std::to_string(mpi::size())+".py");
    if (mpi::size() == 4) {
        const auto& polygons = input_fs.polygons();
        double tol_xy = 0.05; // in metres!
        EXPECT_APPROX_EQ(polygons[0].xy(),(std::vector<PointXY>{{-60000,-85003.7},{1250.01,-85003.7},{1250.01,-71253.7},{-1249.99,-71253.7},{-1249.99,1246.27},{-60000,1246.27},{-60000,-85003.7}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[1].xy(),(std::vector<PointXY>{{1250.01,-85003.7},{60000,-85003.7},{60000,-1253.73},{1250.01,-1253.73},{1250.01,1246.27},{-1249.99,1246.27},{-1249.99,-71253.7},{1250.01,-71253.7},{1250.01,-85003.7}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[2].xy(),(std::vector<PointXY>{{-60000,1246.27},{1250.01,1246.27},{1250.01,73746.3},{-1249.99,73746.3},{-1249.99,84996.3},{-60000,84996.3},{-60000,1246.27}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[3].xy(),(std::vector<PointXY>{{1250.01,-1253.73},{60000,-1253.73},{60000,84996.3},{-1249.99,84996.3},{-1249.99,73746.3},{1250.01,73746.3},{1250.01,-1253.73}}), tol_xy);
    }
}

void check_output_structured_partition_polygon(FunctionSpace output_fs, const std::string& name) {
    ATLAS_TRACE();
    output_fs.polygon().outputPythonScript("output_partitions_n"+std::to_string(mpi::size())+"_"+name+".py");
    if (name == "lonlat_centred" && mpi::size() == 4) {
        const auto& polygons = output_fs.polygons();
        double tol_xy = 0.0005; // in degrees!
        EXPECT_APPROX_EQ(polygons[0].xy(),(std::vector<PointXY>{{4.4,51.467},{4.8875,51.467},{4.8875,51.9795},{4.4,51.9795},{4.4,51.467}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[1].xy(),(std::vector<PointXY>{{4.8875,51.467},{5.4,51.467},{5.4,51.9545},{4.9125,51.9545},{4.9125,51.9795},{4.8875,51.9795},{4.8875,51.467}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[2].xy(),(std::vector<PointXY>{{4.4,51.9795},{4.9125,51.9795},{4.9125,52.467},{4.4,52.467},{4.4,51.9795}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[3].xy(),(std::vector<PointXY>{{4.9125,51.9545},{5.4,51.9545},{5.4,52.467},{4.9125,52.467},{4.9125,51.9545}}), tol_xy);
    }
    else if (name == "lonlat_bottomleft" && mpi::size() == 4) {
        const auto& polygons = output_fs.polygons();
        double tol_xy = 0.0005; // in degrees!
        EXPECT_APPROX_EQ(polygons[0].xy(),(std::vector<PointXY>{{4.2,51.3},{4.7,51.3},{4.7,51.8},{4.2,51.8},{4.2,51.3}}), tol_xy);
        EXPECT_APPROX_EQ(polygons[1].xy(),(std::vector<PointXY>{}), tol_xy);
        EXPECT_APPROX_EQ(polygons[2].xy(),(std::vector<PointXY>{}), tol_xy);
        EXPECT_APPROX_EQ(polygons[3].xy(),(std::vector<PointXY>{}), tol_xy);
    }
}

void check_input_mesh_partition_polygon(FunctionSpace input_fs) {
    ATLAS_TRACE();
    input_fs.polygon().outputPythonScript("input_mesh_partitions_n"+std::to_string(mpi::size())+".py");
    if (mpi::size() == 4) {
        [[maybe_unused]] const auto& polygons = input_fs.polygons();
        // double tol_xy = 0.05; // in metres!
        // EXPECT_APPROX_EQ(polygons[0].xy(),(std::vector<PointXY>{{-60000,-85003.7},{1250.01,-85003.7},{1250.01,-71253.7},{-1249.99,-71253.7},{-1249.99,1246.27},{-60000,1246.27},{-60000,-85003.7}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[1].xy(),(std::vector<PointXY>{{1250.01,-85003.7},{60000,-85003.7},{60000,-1253.73},{1250.01,-1253.73},{1250.01,1246.27},{-1249.99,1246.27},{-1249.99,-71253.7},{1250.01,-71253.7},{1250.01,-85003.7}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[2].xy(),(std::vector<PointXY>{{-60000,1246.27},{1250.01,1246.27},{1250.01,73746.3},{-1249.99,73746.3},{-1249.99,84996.3},{-60000,84996.3},{-60000,1246.27}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[3].xy(),(std::vector<PointXY>{{1250.01,-1253.73},{60000,-1253.73},{60000,84996.3},{-1249.99,84996.3},{-1249.99,73746.3},{1250.01,73746.3},{1250.01,-1253.73}}), tol_xy);
    }
}

void check_output_mesh_partition_polygon(FunctionSpace output_fs, const std::string& name) {
    ATLAS_TRACE();
    output_fs.polygon().outputPythonScript("output_mesh_partitions_n"+std::to_string(mpi::size())+"_"+name+".py");
    if (name == "lonlat_centred" && mpi::size() == 4) {
        [[maybe_unused]] const auto& polygons = output_fs.polygons();
        // double tol_xy = 0.0005; // in degrees!
        // EXPECT_APPROX_EQ(polygons[0].xy(),(std::vector<PointXY>{{4.4,51.467},{4.8875,51.467},{4.8875,51.9795},{4.4,51.9795},{4.4,51.467}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[1].xy(),(std::vector<PointXY>{{4.8875,51.467},{5.4,51.467},{5.4,51.9545},{4.9125,51.9545},{4.9125,51.9795},{4.8875,51.9795},{4.8875,51.467}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[2].xy(),(std::vector<PointXY>{{4.4,51.9795},{4.9125,51.9795},{4.9125,52.467},{4.4,52.467},{4.4,51.9795}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[3].xy(),(std::vector<PointXY>{{4.9125,51.9545},{5.4,51.9545},{5.4,52.467},{4.9125,52.467},{4.9125,51.9545}}), tol_xy);
    }
    else if (name == "lonlat_bottomleft" && mpi::size() == 4) {
        [[maybe_unused]] const auto& polygons = output_fs.polygons();
        // double tol_xy = 0.0005; // in degrees!
        // EXPECT_APPROX_EQ(polygons[0].xy(),(std::vector<PointXY>{{4.2,51.3},{4.7,51.3},{4.7,51.8},{4.2,51.8},{4.2,51.3}}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[1].xy(),(std::vector<PointXY>{}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[2].xy(),(std::vector<PointXY>{}), tol_xy);
        // EXPECT_APPROX_EQ(polygons[3].xy(),(std::vector<PointXY>{}), tol_xy);
    }
}



FunctionSpace output_functionspace( const std::string& name, const FunctionSpace& input_functionspace) {
    ATLAS_TRACE();
    Grid grid(get_output_grid(name));
    auto partitioner = grid::MatchingPartitioner(input_functionspace);
    grid::Distribution dist = partitioner.partition(grid);
    return functionspace::StructuredColumns(grid, dist);
}


FieldSet create_source_fields(FunctionSpace& fs, idx_t nb_fields, idx_t nb_levels) {
    using Value = double;
    FieldSet fields_source;
    auto lonlat = array::make_view<double, 2>(fs.lonlat());
    for (idx_t f = 0; f < nb_fields; ++f) {
        auto field_source = fields_source.add(fs.createField<Value>());
        auto source       = array::make_view<Value, 2>(field_source);
        for (idx_t n = 0; n < fs.size(); ++n) {
            for (idx_t k = 0; k < nb_levels; ++k) {
                source(n, k) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 100. + double(k) / 2);
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

void do_structured_test( std::string type, int input_halo, bool matrix_free, const std::string& output_grid_name ) {
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid(get_input_grid());
    output::Gmsh gmsh("input_mesh.msh", util::Config("coordinates","lonlat"));
    gmsh.write(Mesh(input_grid));

    {
        output::Gmsh gmsh("input_mesh_xy.msh", util::Config("coordinates","xy"));
        gmsh.write(Mesh(input_grid));
    }

    functionspace::StructuredColumns input_fs(input_grid, option::levels(nb_levels) |
        option::halo(input_halo));

    check_input_structured_partition_polygon(input_fs);

    FieldSet fields_source = create_source_fields(input_fs, nb_fields, nb_levels);
    gmsh.write(fields_source);

    FunctionSpace output_fs = output_functionspace(output_grid_name, input_fs);

    check_output_structured_partition_polygon(output_fs, output_grid_name);

    output::Gmsh output_gmsh("output_mesh.msh", util::Config("coordinates","lonlat"));
    output_gmsh.write( MeshGenerator("structured").generate(get_output_grid(output_grid_name),grid::MatchingPartitioner(input_fs)));

    Interpolation interpolation(option::type(type) |
        util::Config("matrix_free",matrix_free) |
        util::Config("verbose",eckit::Resource<bool>("--verbose",false)),
        input_fs, output_fs);

    FieldSet fields_target = create_target_fields(output_fs, nb_fields, nb_levels);

    interpolation.execute(fields_source, fields_target);

    output_gmsh.write(fields_target);
}

void do_nodes_test( std::string type, int input_halo, bool matrix_free, const std::string& output_grid_name ) {
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid(get_input_grid());
    Mesh input_mesh(input_grid);
    output::Gmsh gmsh("input_mesh.msh", util::Config("coordinates","lonlat"));
    gmsh.write(input_mesh);

    {
        output::Gmsh gmsh("input_mesh_xy.msh", util::Config("coordinates","xy"));
        gmsh.write(input_mesh);
    }

    functionspace::NodeColumns input_fs(input_mesh, option::levels(nb_levels) |
        option::halo(input_halo));

    check_input_mesh_partition_polygon(input_fs);

    FieldSet fields_source = create_source_fields(input_fs, nb_fields, nb_levels);
    gmsh.write(fields_source);

    FunctionSpace output_fs = output_functionspace(output_grid_name, input_fs);

    check_output_mesh_partition_polygon(output_fs, output_grid_name);

    output::Gmsh output_gmsh("output_mesh_"+output_grid_name+"_"+type+".msh", util::Config("coordinates","lonlat"));
    output_gmsh.write( MeshGenerator("structured").generate(get_output_grid(output_grid_name),grid::MatchingPartitioner(input_fs)));

    Interpolation interpolation(option::type(type) |
        util::Config("matrix_free",matrix_free) |
        util::Config("verbose",eckit::Resource<bool>("--verbose",false)),
        input_fs, output_fs);

    FieldSet fields_target = create_target_fields(output_fs, nb_fields, nb_levels);

    interpolation.execute(fields_source, fields_target);

    output_gmsh.write(fields_target);
}

// For now regional-bilinear and (k-)nearest-neighbour is the only method that works with StructuredColumns.
// There are still issues with StructuredInterpolation methods (cubic, ...) which expect a global grid
CASE("test structured regional-bilinear, output_grid=lonlat_centred, halo 0") {
    EXPECT_NO_THROW( do_structured_test("regional-bilinear",0,false, "lonlat_centred") );
}

CASE("test structured regional-bilinear, output_grid=lonlat_bottomleft, halo 0") {
    // With MPI_SIZE=4, only rank 0 will have output grid points in this case, which is also a good test
    EXPECT_NO_THROW( do_structured_test("regional-bilinear",0,false, "lonlat_bottomleft") );
}

CASE("test structured nearest-neighbour, output_grid=lonlat_centred, halo 0") {
    EXPECT_NO_THROW( do_structured_test("nearest-neighbour",0,false, "lonlat_centred") );
}

CASE("test structured nearest-neighbour, output_grid=lonlat_bottomleft, halo 0") {
    // With MPI_SIZE=4, only rank 0 will have output grid points in this case, which is also a good test
    EXPECT_NO_THROW( do_structured_test("nearest-neighbour",0,false, "lonlat_bottomleft") );
}

// Interpolation methods that require the mesh
CASE("test meshed finite-element, output_grid=lonlat_centred, halo 0") {
    EXPECT_NO_THROW( do_nodes_test("finite-element",0,false, "lonlat_centred") );
}

CASE("test meshed finite-element, output_grid=lonlat_bottomleft, halo 0") {
    // With MPI_SIZE=4, only rank 0 will have output grid points in this case, which is also a good test
    EXPECT_NO_THROW( do_nodes_test("finite-element",0,false, "lonlat_bottomleft") );
}

CASE("test meshed nearest-neighbour, output_grid=lonlat_centred, halo 0") {
    EXPECT_NO_THROW( do_nodes_test("nearest-neighbour",0,false, "lonlat_centred") );
}

CASE("test meshed nearest-neighbour, output_grid=lonlat_bottomleft, halo 0") {
    // With MPI_SIZE=4, only rank 0 will have output grid points in this case, which is also a good test
    EXPECT_NO_THROW( do_nodes_test("nearest-neighbour",0,false, "lonlat_bottomleft") );
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
