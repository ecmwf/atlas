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
#include <fstream>
#include <string>

#include "eckit/types/FloatCompare.h"

#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/grid/Unstructured.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/detail/PointCloudIO.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace {

using mdspan_xy = atlas::mdspan<const double, atlas::extents<size_t, atlas::dynamic_extent, 2>>;
mdspan_xy make_mdspan(const atlas::Field& xy) {
    return mdspan_xy{xy.array().host_data<const double>(), xy.shape(0), 2 };
}

}

namespace test_arrays {

const size_t nb_pts     = 5;
const size_t nb_fld     = 2;
const size_t nb_columns = 2 + nb_fld;

const double lon[] = {-31.233, -28.717, -27.217, -25.750, -16.917};
const double lat[] = {39.467, 38.583, 38.483, 37.817, 32.650};

const double helper_f1[] = {1., 2., 3., 4., 5.};
const double helper_f2[] = {-0.1, -0.2, -0.3, -0.4, -0.5};

const double* fvalues[] = {helper_f1, helper_f2};
const char* fnames[]    = {" f_1  ", "f    2 "};

}  // namespace test_arrays

namespace test_vectors {

const size_t nb_pts     = test_arrays::nb_pts;
const size_t nb_fld     = test_arrays::nb_fld;
const size_t nb_columns = test_arrays::nb_columns;

const std::vector<double> lon(test_arrays::lon, test_arrays::lon + nb_pts);
const std::vector<double> lat(test_arrays::lat, test_arrays::lat + nb_pts);

std::vector<double> helper_f1(test_arrays::helper_f1, test_arrays::helper_f1 + nb_pts);
std::vector<double> helper_f2(test_arrays::helper_f2, test_arrays::helper_f2 + nb_pts);
std::vector<double>* helper_fvalues[] = {&helper_f1, &helper_f2};

const std::vector<std::vector<double>*> fvalues(helper_fvalues, helper_fvalues + nb_fld);
const std::vector<std::string> fnames(test_arrays::fnames, test_arrays::fnames + nb_fld);

}  // namespace test_vectors

void test_write_file(const std::string& file_path, const size_t& nb_pts, const size_t& nb_columns) {
    REQUIRE(nb_pts > 0);
    std::ofstream f(file_path.c_str());
    if (!f.is_open()) {
        throw eckit::CantOpenFile(file_path);
    }
    f << "PointCloudIO " << nb_pts << "	" << nb_columns
      << "  lon	lat	f_1				"
         "__f3			more resilience	\n"
         "-31.233	39.467	1.	-0.1\n"
         "-28.717	38.583	2.	-0.2 even	more "
         "resilience\n"
         "-27.217	38.483	3.	-0.3\n"
         "-25.750	37.817	4.	-0.4\n"
         "-16.917	32.650	5.	-0.5\n";
}

void test_write_file_bad(const std::string& file_path) {
    std::ofstream f(file_path.c_str());
    if (!f.is_open()) {
        throw eckit::CantOpenFile(file_path);
    }
    f << '?';
}

}  // end anonymous namespace

// ------------------------------------------------------------------

namespace atlas {
namespace test {

// ------------------------------------------------------------------

CASE("read_inexistent_file") {
    EXPECT_THROWS_AS(output::detail::PointCloudIO::read("pointcloud.txt_should_not_exist"), eckit::CantOpenFile);
}

CASE("read_badly_formatted_file") {
    test_write_file_bad("pointcloud.txt");
    EXPECT_THROWS_AS(output::detail::PointCloudIO::read("pointcloud.txt"), eckit::Exception);
}

CASE("read_grid_sample_file") {
    // test sample file, header properly formatted (some fluff is present)
    test_write_file("pointcloud.txt", test_arrays::nb_pts, test_arrays::nb_columns);

    Log::info() << "pointcloud.txt created" << std::endl;

    Mesh mesh = output::detail::PointCloudIO::read("pointcloud.txt");

    Log::info() << "Mesh created" << std::endl;
    Grid grid(new grid::detail::grid::Unstructured(make_mdspan(mesh.nodes().xy())));
    EXPECT(grid);

    EXPECT(grid.size() == test_arrays::nb_pts);
    EXPECT(mesh.nodes().has_field("f_1") == true);
    EXPECT(mesh.nodes().has_field("f3") == true);
}

CASE("read_grid_sample_file_header_less_rows") {
    Log::info() << "Creating Mesh..." << std::endl;
    // test sample file with (wrong) header with less rows
    test_write_file("pointcloud.txt", test_arrays::nb_pts - 2, test_arrays::nb_columns);

    Log::info() << "Creating Mesh..." << std::endl;
    Mesh mesh = output::detail::PointCloudIO::read("pointcloud.txt");
    Log::info() << "Creating Mesh...done" << std::endl;
    Grid grid(new grid::detail::grid::Unstructured(make_mdspan(mesh.nodes().xy())));
    EXPECT(grid);

    EXPECT(grid.size() == test_arrays::nb_pts - 2);
    EXPECT(mesh.nodes().has_field("f_1") == true);
    EXPECT(mesh.nodes().has_field("f3") == true);
}

CASE("read_grid_sample_file_header_less_columns_1") {
    // test sample file with (wrong) header with one field less
    test_write_file("pointcloud.txt", test_arrays::nb_pts, test_arrays::nb_columns - 1);

    Mesh mesh = output::detail::PointCloudIO::read("pointcloud.txt");
    Grid grid(new grid::detail::grid::Unstructured(make_mdspan(mesh.nodes().xy())));
    EXPECT(grid);

    EXPECT(grid.size() == test_arrays::nb_pts);
    EXPECT(mesh.nodes().has_field("f_1") == true);
    EXPECT(mesh.nodes().has_field("f3") == false);
}

CASE("read_grid_sample_file_header_less_columns_2") {
    // test sample file with (wrong) header with no fields
    test_write_file("pointcloud.txt", test_arrays::nb_pts, test_arrays::nb_columns - test_arrays::nb_fld);

    Mesh mesh = output::detail::PointCloudIO::read("pointcloud.txt");
    Grid grid(new grid::detail::grid::Unstructured(make_mdspan(mesh.nodes().xy())));
    EXPECT(grid);

    EXPECT(grid.size() == test_arrays::nb_pts);
    EXPECT(mesh.nodes().has_field("f_1") == false);
    EXPECT(mesh.nodes().has_field("f3") == false);
}

CASE("write_array") {
    std::ifstream f;
    std::string signature, str_lon, str_lat, str_f1, str_f2;
    size_t nb_pts, nb_columns;

    output::detail::PointCloudIO::write("pointcloud.txt", test_arrays::nb_pts, test_arrays::lon, test_arrays::lat);
    f.open("pointcloud.txt");
    EXPECT(f.is_open());
    f >> signature >> nb_pts >> nb_columns >> str_lon >> str_lat >> str_f1 >> str_f2;
    f.close();

    EXPECT(nb_pts == test_arrays::nb_pts);
    EXPECT(nb_columns == test_arrays::nb_columns - test_arrays::nb_fld);
    EXPECT(str_lon == "lon");
    EXPECT(str_lat == "lat");
    EXPECT(str_f1 != "f_1");     // (this column is not written)
    EXPECT(str_f2 != "f____2");  // (this column is not written)
}

CASE("write_array_less_rows") {
    std::ifstream f;
    std::string signature, str_lon, str_lat, str_f1, str_f2;
    size_t nb_pts, nb_columns;

    output::detail::PointCloudIO::write("pointcloud.txt", test_arrays::nb_pts - 1 /* deliberate */, test_arrays::lon,
                                        test_arrays::lat, test_arrays::nb_fld, test_arrays::fvalues,
                                        test_arrays::fnames);
    f.open("pointcloud.txt");
    EXPECT(f.is_open());
    f >> signature >> nb_pts >> nb_columns >> str_lon >> str_lat >> str_f1 >> str_f2;
    f.close();

    EXPECT(nb_pts == test_arrays::nb_pts - 1);  // (one row is not written)
    EXPECT(nb_columns == test_arrays::nb_columns);
    EXPECT(str_lon == "lon");
    EXPECT(str_lat == "lat");
    EXPECT(str_f1 == "f_1");
    EXPECT(str_f2 == "f____2");
}

CASE("write_array_less_columns") {
    std::ifstream f;
    std::string signature, str_lon, str_lat, str_f1, str_f2;
    size_t nb_pts, nb_columns;

    output::detail::PointCloudIO::write("pointcloud.txt", test_arrays::nb_pts, test_arrays::lon, test_arrays::lat,
                                        test_arrays::nb_fld - 1 /* deliberate */, test_arrays::fvalues,
                                        test_arrays::fnames);
    f.open("pointcloud.txt");
    EXPECT(f.is_open());
    f >> signature >> nb_pts >> nb_columns >> str_lon >> str_lat >> str_f1 >> str_f2;
    f.close();

    EXPECT(nb_pts == test_arrays::nb_pts);
    EXPECT(nb_columns == test_arrays::nb_columns - 1);  // (one column is not written)
    EXPECT(str_lon == "lon");
    EXPECT(str_lat == "lat");
    EXPECT(str_f1 == "f_1");
    EXPECT(str_f2 != "f____2");  // (this column is not written)
}

CASE("write_vector_all_fields") {
    std::ifstream f;
    std::string signature, str_lon, str_lat, str_f1, str_f2;
    size_t nb_pts, nb_columns;

    output::detail::PointCloudIO::write("pointcloud.txt", test_vectors::lon, test_vectors::lat, test_vectors::fvalues,
                                        test_vectors::fnames);
    f.open("pointcloud.txt");
    EXPECT(f.is_open());
    f >> signature >> nb_pts >> nb_columns >> str_lon >> str_lat >> str_f1 >> str_f2;
    f.close();

    EXPECT(nb_pts == test_vectors::nb_pts);
    EXPECT(nb_columns == test_vectors::nb_columns);
    EXPECT(str_lon == "lon");
    EXPECT(str_lat == "lat");
    EXPECT(str_f1 == "f_1");
    EXPECT(str_f2 == "f____2");
}

CASE("write_vector_no_fields") {
    std::ifstream f;
    std::string signature, str_lon, str_lat, str_f1, str_f2;
    size_t nb_pts, nb_columns;

    output::detail::PointCloudIO::write("pointcloud.txt", test_vectors::lon, test_vectors::lat);
    f.open("pointcloud.txt");
    EXPECT(f.is_open());
    f >> signature >> nb_pts >> nb_columns >> str_lon >> str_lat >> str_f1 >> str_f2;
    f.close();

    EXPECT(nb_pts == test_vectors::nb_pts);
    EXPECT(nb_columns == test_vectors::nb_columns - test_vectors::nb_fld);
    EXPECT(str_lon == "lon");
    EXPECT(str_lat == "lat");
    EXPECT(str_f1 != "f_1");     // (this column is not written)
    EXPECT(str_f2 != "f____2");  // (this column is not written)
}

static double funny_formula(int x) {
    return ((double)(x)) * std::pow((double)-1., (int)(x));
}

CASE("write_read_write_field") {
    // build suitable data structures do hold field name & values
    std::string field_name("my_super_field");
    std::vector<double> field_values(test_vectors::nb_pts, 0.);
    for (size_t i = 0; i < test_vectors::nb_pts; ++i) {
        field_values[i] = funny_formula(i);
    }

    // PART 1
    // write field vector values as column in file "pointcloud.txt"
    Log::info() << "Part 1" << std::endl;

    std::ifstream f;
    std::string signature, str_lon, str_lat, str_f;
    size_t nb_pts, nb_columns;

    output::detail::PointCloudIO::write("pointcloud.txt", test_vectors::lon, test_vectors::lat,
                                        std::vector<std::vector<double>*>(1, &field_values),
                                        std::vector<std::string>(1, field_name));
    f.open("pointcloud.txt");
    EXPECT(f.is_open());
    f >> signature >> nb_pts >> nb_columns >> str_lon >> str_lat >> str_f;
    f.close();

    EXPECT(nb_pts == test_vectors::nb_pts);
    EXPECT(nb_columns == 2 + 1);  // (lon,lat,my_super_field)
    EXPECT(str_lon == "lon");
    EXPECT(str_lat == "lat");
    EXPECT(str_f == "my_super_field");

    // PART 2
    // read field vector from just-created file
    Log::info() << "Part 2" << std::endl;

    Mesh mesh = output::detail::PointCloudIO::read("pointcloud.txt");
    Grid grid(new grid::detail::grid::Unstructured(make_mdspan(mesh.nodes().xy())));
    EXPECT(grid);

    EXPECT(grid.size() == test_vectors::nb_pts);

    mesh::Nodes& nodes = mesh.nodes();
    EXPECT(nodes.has_field("my_super_field") == true);
    EXPECT(nodes.has_field("_StRaNgE_FiElD_NaMe_") == false);

    // PART 3
    // check field values to a very small tolerance (relative tol. 0.001%)
    Log::info() << "Part 3" << std::endl;

    Field& field(nodes.field("my_super_field"));
    EXPECT(
        /* data used to write file*/ test_vectors::nb_pts ==
        /* data read from file*/ field.size());

    array::ArrayView<double, 1> field_data = array::make_view<double, 1>(field);
    for (idx_t i = 0; i < field_data.size(); ++i) {
        EXPECT(eckit::types::is_approximately_equal(funny_formula(i), field_data(i),
                                                    0.001));  // 0.001% relative error
        EXPECT(eckit::types::is_approximately_equal(funny_formula(i), field_data(i),
                                                    0.001));  // 0.001% relative error
    }

    // PART 4
    // write to file a Field (the just-read one),
    // a FieldSet, and
    // a Grid (should be exactly the same)
    Log::info() << "Part 4" << std::endl;

    FieldSet fieldset;
    EXPECT_NO_THROW(fieldset.add(field));

    functionspace::NodeColumns functionspace(mesh);

    EXPECT_NO_THROW(output::detail::PointCloudIO::write("pointcloud_FieldSet.txt", fieldset, functionspace));
    EXPECT_NO_THROW(output::detail::PointCloudIO::write("pointcloud_Grid.txt", mesh));

    Mesh mesh_from_FieldSet = output::detail::PointCloudIO::read("pointcloud_FieldSet.txt");
    Grid grid_from_FieldSet(new grid::detail::grid::Unstructured(make_mdspan(mesh_from_FieldSet.nodes().xy())));

    Mesh mesh_from_Grid(output::detail::PointCloudIO::read("pointcloud_Grid.txt"));
    Grid grid_from_Grid(new grid::detail::grid::Unstructured(make_mdspan(mesh_from_Grid.nodes().xy())));

    EXPECT(grid_from_FieldSet);
    EXPECT(grid_from_Grid);


    // PART 5
    // compare reading of reference data to:
    // - grid_from_FieldSet, and
    // - grid_from_Grid (all different but equivalent writing methods)
    Log::info() << "Part 5" << std::endl;

    // (header section)
    EXPECT(grid_from_FieldSet.size() == test_arrays::nb_pts);
    EXPECT(mesh_from_FieldSet.nodes().has_field("my_super_field") == true);
    EXPECT(mesh_from_FieldSet.nodes().has_field("_StRaNgE_FiElD_NaMe_") == false);

    EXPECT(grid_from_Grid.size() == test_arrays::nb_pts);
    EXPECT(mesh_from_FieldSet.nodes().has_field("my_super_field") == true);
    EXPECT(mesh_from_FieldSet.nodes().has_field("_StRaNgE_FiElD_NaMe_") == false);

    // (data section: guarantee data are from different places, to make checks
    // useful)
    const Field& field_from_FieldSet(mesh_from_FieldSet.nodes().field("my_super_field"));
    const Field& field_from_Grid(mesh_from_FieldSet.nodes().field("my_super_field"));
    EXPECT(field.array().data<double>() != field_from_FieldSet.array().data<double>());
    EXPECT(field.array().data<double>() != field_from_Grid.array().data<double>());

    auto field_from_FieldSet_data = array::make_view<double, 1>(field_from_FieldSet);
    auto field_from_Grid_data     = array::make_view<double, 1>(field_from_Grid);
    for (size_t i = 0; i < test_arrays::nb_pts; ++i) {
        EXPECT(eckit::types::is_approximately_equal(field_data(i), field_from_FieldSet_data(i),
                                                    0.001));  // 0.001% relative error
        EXPECT(eckit::types::is_approximately_equal(field_data(i), field_from_Grid_data(i), 0.001));  // ...
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
