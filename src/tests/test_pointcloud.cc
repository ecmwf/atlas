/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestPointCloud
#include "ecbuild/boost_test_framework.h"


#include <string>

#include "atlas/atlas_config.h"

#include "atlas/grids/Unstructured.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/io/PointCloud.h"


namespace {


  bool pointcloud_write_test_file(const std::string& file_path, int nb_pts=5, int nb_columns=4)
  {
    if (nb_pts<1 || nb_columns<2)
      return false;
    std::ofstream f(file_path.c_str());
    return (f?
              f << "PointCloud " << nb_pts << "	" << nb_columns << "  lon	lat	f_1				__f3			more resilience	\n"
                   "-31.233	39.467	1.	-0.1\n"
                   "-28.717	38.583	2.	-0.2 even	more resilience\n"
                   "-27.217	38.483	3.	-0.3\n"
                   "-25.750	37.817	4.	-0.4\n"
                   "-16.917	32.650	5.	-0.5\n"
            : false);
  }


  namespace pointcloud_test_arrays {


    const double
        lon[] = { -31.233, -28.717, -27.217, -25.750, -16.917 },
        lat[] = {  39.467,  38.583,  38.483,  37.817,  32.650 },
        f1[] = {  1.,   2.,   3.,   4.,   5.  },
        f2[] = { -0.1, -0.2, -0.3, -0.4, -0.5 };
    const double *afvalues[] = {  f1,   f2  };
    const char   *afnames [] = { " f_1  ", "f    2 " };


  }
}


using namespace atlas;


BOOST_AUTO_TEST_CASE( test_pointcloud_read_grid_sample_file )
{
  // test sample file, header properly formatted (some fluff is present)
  BOOST_REQUIRE(pointcloud_write_test_file("pointcloud.txt"));

  grids::Unstructured* grid = io::PointCloud::read("pointcloud.txt");
  BOOST_REQUIRE(grid);

  BOOST_CHECK_EQUAL(grid->npts(),5);
  BOOST_CHECK_EQUAL(true,grid->mesh().has_function_space("nodes"));
  BOOST_CHECK_EQUAL(true,grid->mesh().function_space(0).has_field("f_1"));
  BOOST_CHECK_EQUAL(true,grid->mesh().function_space(0).has_field("f3"));

  delete grid;
}


BOOST_AUTO_TEST_CASE( test_pointcloud_read_grid_sample_file_header_less_rows )
{
  // test sample file with (wrong) header with less rows
  BOOST_REQUIRE(pointcloud_write_test_file("pointcloud.txt",3,4));

  grids::Unstructured* grid = io::PointCloud::read("pointcloud.txt");
  BOOST_REQUIRE(grid);

  BOOST_CHECK_EQUAL(grid->npts(),3);
  BOOST_CHECK_EQUAL(true,grid->mesh().function_space(0).has_field("f_1"));
  BOOST_CHECK_EQUAL(true,grid->mesh().function_space(0).has_field("f3"));

  delete grid;
}


BOOST_AUTO_TEST_CASE( test_pointcloud_read_grid_sample_file_header_less_columns )
{
  // test sample file with (wrong) header with less columns
  BOOST_REQUIRE(pointcloud_write_test_file("pointcloud.txt",5,3));

  grids::Unstructured* grid = io::PointCloud::read("pointcloud.txt");
  BOOST_REQUIRE(grid);

  BOOST_CHECK_EQUAL(grid->npts(),5);
  BOOST_CHECK_EQUAL(true, grid->mesh().function_space(0).has_field("f_1"));
  BOOST_CHECK_EQUAL(false,grid->mesh().function_space(0).has_field("f3"));

  delete grid;
  BOOST_REQUIRE(pointcloud_write_test_file("pointcloud.txt",5,2));

  grid = io::PointCloud::read("pointcloud.txt");
  BOOST_REQUIRE(grid);

  BOOST_CHECK_EQUAL(grid->npts(),5);
  BOOST_CHECK_EQUAL(false,grid->mesh().function_space(0).has_field("f_1"));
  BOOST_CHECK_EQUAL(false,grid->mesh().function_space(0).has_field("f3"));

  delete grid;
}


BOOST_AUTO_TEST_CASE( test_pointcloud_write_array )
{
  using namespace atlas;

  using pointcloud_test_arrays::lon;
  using pointcloud_test_arrays::lat;
  using pointcloud_test_arrays::afvalues;
  using pointcloud_test_arrays::afnames;

  std::ifstream f;
  std::string signature;
  size_t
      nb_pts,
      nb_columns;

  io::PointCloud::write("pointcloud_0f.txt", 5, lon, lat);
  f.open("pointcloud_0f.txt");
  f >> signature >> nb_pts >> nb_columns;
  f.close();
  BOOST_CHECK_EQUAL(5,nb_pts);
  BOOST_CHECK_EQUAL(2,nb_columns);

  io::PointCloud::write("pointcloud_1f.txt", 5, lon, lat, 1, afvalues, afnames);
  f.open("pointcloud_1f.txt");
  f >> signature >> nb_pts >> nb_columns;
  f.close();
  BOOST_CHECK_EQUAL(5,nb_pts);
  BOOST_CHECK_EQUAL(3,nb_columns);

  io::PointCloud::write("pointcloud_2f.txt", 4, lon, lat, 2, afvalues, afnames);
  f.open("pointcloud_2f.txt");
  f >> signature >> nb_pts >> nb_columns;
  f.close();
  BOOST_CHECK_EQUAL(4,nb_pts);
  BOOST_CHECK_EQUAL(4,nb_columns);
}


