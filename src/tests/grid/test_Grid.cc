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

#define BOOST_TEST_MODULE TestGrid
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "eckit/io/StdFile.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/GribGridBuilder.h"
#include "atlas/grid/StackGribFile.h"
#include "atlas/grid/GribWrite.h"
#include "atlas/grid/GridSpec.h"


using namespace std;
using namespace eckit;
using namespace atlas::grid;
using namespace atlas;

/// Test for Grid* derivatives
/// This test uses the grib samples directory.
/// We open each file and attempt to create the Grid* class.
/// However this has exposed errors in the grib samples files,
/// especially for reduced gaussian grids. Hence we will need to wait till Shahram has
/// has a chance to fix the issues, before all the tests pass.

static void test_grids_from_grib_sample_directory( const std::string& directory);
static void test_grib_file(const std::string& file);

BOOST_AUTO_TEST_SUITE( TestGrid )

BOOST_AUTO_TEST_CASE( test_grids_from_samples_dir )
{
   cout << "Grid:: ...test_grids_from_samples_dir\n";
   BOOST_CHECK(true); // stop boost test from complaining about no checks

   // Traverse all the GRIB samples files, for gridType first determine sample dir
   std::vector<std::string> sample_dirs;
   GribWrite::determine_grib_samples_dir(sample_dirs);
   BOOST_REQUIRE_MESSAGE(!sample_dirs.empty(),"Expected sample dirs to be found");

   // now test these dirs
   for(size_t i = 0; i < sample_dirs.size(); ++i) {
      test_grids_from_grib_sample_directory(sample_dirs[i]);
   }
}

BOOST_AUTO_TEST_SUITE_END()


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
   StackGribFile the_grib_file(the_file_path);

   std::cout << " Get the grid type" << std::endl;
   char string_value[64];
   size_t len = sizeof(string_value)/sizeof(char);
   int err = grib_get_string(the_grib_file.handle(),"gridType",string_value,&len);
   BOOST_CHECK_MESSAGE(err == 0,"grib_get_string(gridType) failed for \nfile " << the_file_path << " IGNORING !!!!\n");

   std::string gridType = string_value;
   std::cout << " Create Grid derivatives " << gridType << std::endl;
   if ( gridType == "polar_stereographic" || gridType == "sh" ) {
      std::cout << " ** Ignoring grid types [ polar_stereographic | sh ] " << std::endl;
      return;
   }

   long editionNumber = 0;
   GRIB_CHECK(grib_get_long(the_grib_file.handle(),"editionNumber",&editionNumber),0);


   // Unstructured grid can not handle Spherical harmonics
   atlas::grid::Grid::Ptr the_grid = GRIBGridBuilder::instance().build_grid_from_grib_handle(the_grib_file.handle());
   BOOST_CHECK_MESSAGE(the_grid,"GRIBGridBuilder::instance().build_grid_from_grib_handle failed for file " << the_file_path);
   if (!the_grid) return;

   // The Grid produced, has a GRID spec, the grid spec can be used to,
   // make sure the grid types match
   eckit::ScopedPtr< GridSpec > the_grid_spec( the_grid->spec() );
   // std::cout << *the_grid_spec << "\n";  // this will kill the performance since points are also written

   BOOST_CHECK_MESSAGE(the_grid->gridType() == gridType,"gridType(" << gridType << ") did not match Grid constructor(" << the_grid->gridType() << ") for file " << the_file_path);
   BOOST_CHECK_MESSAGE(the_grid_spec->grid_type() == gridType,"gridType(" << gridType << ") did not match GridSpec constructor(" << the_grid_spec->grid_type() << ") for file " << the_file_path);

   // find the corresponding sample file.
   // However we need to take into account that the GRIB samples, file are *NOT* unique in their GRID definition.
   // The sample file name produced does not have '.tmpl' extension
   std::string generated_sample_file_name = GribWrite::grib_sample_file( *the_grid_spec , editionNumber);
   BOOST_CHECK_MESSAGE( !generated_sample_file_name.empty()," Could *not* find sample file for grid_spec " << *the_grid_spec );


   // Note: many of the grib samples files are not UNIQUE in their grid specification:
   // hence the use of WARN.
   // remove .tmpl and get base part
   eckit::LocalPathName path(the_file_path);
   LocalPathName the_base_name = path.baseName(false);
   std::string grib_sample_file = the_base_name.localPath();
   BOOST_WARN_MESSAGE( generated_sample_file_name == grib_sample_file, "\nCould not match samples expected '"
                       << grib_sample_file << "' but found('"
                       << generated_sample_file_name
                       << "') for grid spec "
                       << *the_grid_spec );
}
