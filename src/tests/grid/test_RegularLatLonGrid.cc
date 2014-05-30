#define BOOST_TEST_MODULE TestGrid
/*
 * (C) Copyright 1996-2012 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#include <string>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "eckit/io/StdFile.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridBuilder.h"


using namespace std;
using namespace eckit;
using namespace atlas::grid;

/// Test for Grid* derivatives
/// This test uses the grib samples directory.
/// We open each file and attempt to create the Grid* class.
/// However this has exposed errors in the grib samples files,
/// especially for reduced gaussian grids. Hence we will need to wait till Shahram has
/// has a chance to fix the issues, before all the tests pass.

static std::string determine_grib_samples_dir();
static void test_grids_from_grib_sample_directory( const std::string& directory);
static void test_grib_file(const std::string& file);

BOOST_AUTO_TEST_SUITE( TestGrid )

BOOST_AUTO_TEST_CASE( test_grids_from_samples_dir )
{
   cout << "Grid:: ...test_grids_from_samples_dir\n";
   BOOST_CHECK(true); // stop boost test from complaining about no checks

   // Traverse all the GRIB samples files, for gridType first determine sample dir
   std::string samples_dir = determine_grib_samples_dir();
   BOOST_REQUIRE_MESSAGE(!samples_dir.empty(),"Expected sample dirs to be found");

   // now test this dir
   test_grids_from_grib_sample_directory(samples_dir);
}

BOOST_AUTO_TEST_SUITE_END()


static std::string determine_grib_samples_dir()
{
   // Ideally we should use: 'grib_info -t'
   // Try looking for environment variable GRIB_API_INCLUDE
   // GRIB_API_INCLUDE=-I/usr/local/lib/metaps/lib/grib_api/1.10.0/include
   //                  =/usr/local/lib/metaps/lib/grib_api/1.10.0/include /usr/local/apps/jasper/1.900.1/LP64/include /usr/local/apps/jasper/1.900.1/LP64/include
   // samples dir = /usr/local/lib/metaps/lib/grib_api/1.10.0/share/grib_api/samples

   char* include_dir = getenv("GRIB_API_INCLUDE");
   BOOST_REQUIRE_MESSAGE(include_dir,"Expected GRIB_API_INCLUDE to be defined");

   std::string grib_include_dir(include_dir);
   BOOST_REQUIRE_MESSAGE(grib_include_dir.find("grib_api") != std::string::npos,"grib-api not found on directory " << grib_include_dir);

   if (grib_include_dir.find("-I") != std::string::npos) {
      //std::cout << "GRIB_API_INCLUDE=" << grib_include_dir << "\n";
      grib_include_dir.erase(grib_include_dir.begin(),grib_include_dir.begin()+2); // remove -I
   }

   // Handle multiple include dirs
   // If there are any spaces in the string, only take the first include
   size_t space_pos = grib_include_dir.find(" ");
   if (space_pos != std::string::npos) {
      grib_include_dir = grib_include_dir.substr(0,space_pos);
      //std::cout << "GRIB_API_INCLUDE=" << grib_include_dir << "\n";
   }

   // Remove the 'include' and replace with, 'share/grib_api/samples'
   size_t pos = grib_include_dir.find("/include");
   BOOST_REQUIRE_MESSAGE(pos != string::npos,"include not found in directory " << grib_include_dir);

   grib_include_dir = grib_include_dir.replace(pos,grib_include_dir.length(),"/share/grib_api/samples");
   //std::cout << "GRIB_API_INCLUDE=" << grib_include_dir << "\n";

   return grib_include_dir;
}


static void test_grids_from_grib_sample_directory(const std::string& directory)
{
   cout << "*********************************************************************************\n";
   cout << "traversing directory " << directory << "\n";
   PathName dir_path(directory);
   BOOST_CHECK(dir_path.exists());
   BOOST_CHECK(dir_path.isDir());

   int count = 0;
   std::vector<PathName> files;
   std::vector<PathName> directories;
   dir_path.children(files,directories);
   for(size_t i = 0; i < files.size(); i++) {
      try {
         test_grib_file(files[i].localPath());
         //count++;
         //if (count > 2) exit(0);
      }
      catch ( const std::exception & ex )
      {
         std::cout << files[i].localPath() << " " << ex.what() << std::endl;
      }
   }

   // recursively call this function for each directory found
   for(size_t i = 0; i < directories.size(); i++) {
      test_grids_from_grib_sample_directory(directories[i].localPath());
   }
}

static void test_grib_file(const std::string& the_file_path)
{
   std::cout << "\n===================================================================================================" << std::endl;
   std::cout << "Opening GRIB file " << the_file_path << std::endl;
   GribFile the_grib_file(the_file_path);

   std::cout << " Get the grid type" << std::endl;
   char string_value[64];
   size_t len = sizeof(string_value)/sizeof(char);
   int err = grib_get_string(the_grib_file.handle(),"gridType",string_value,&len);
   BOOST_CHECK_MESSAGE(err == 0,"grib_get_string(gridType) failed for \nfile " << the_file_path << " IGNORING !!!!\n");

   std::string gridType = string_value;
   std::cout << " Create Grid derivatives " << gridType << std::endl;

   if ( gridType == "polar_stereographic" || gridType == "rotated_ll" || gridType == "reduced_ll" || gridType == "sh" ) {
      std::cout << " ** Ignoring grid types [ polar_stereographic | rotated_ll | reduced_ll || sh ] " << std::endl;
      return;
   }

   // Unstructured grid can not handle Spherical harmonics
   atlas::grid::Grid::Ptr the_grid = GRIBGridBuilder::instance().build_grid_from_grib_handle(the_grib_file.handle());
   BOOST_CHECK_MESSAGE(the_grid,"GRIBGridBuilder::instance().build_grid_from_grib_handle failed for file " << the_file_path);
   if (!the_grid) return;

   BOOST_CHECK_MESSAGE(the_grid->gridType() == gridType,"gridType(" << gridType << ") did not match Grid constructor(" << the_grid->gridType() << ") for file " << the_file_path);
}
