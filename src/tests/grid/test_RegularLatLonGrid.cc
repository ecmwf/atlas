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

#include "eckit/log/Log.h"
#include <string>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

#include "atlas/grid/Grid.h"
#include "atlas/grid/GribRead.h"
#include "atlas/grid/RegularLatLonGrid.h"
#include "atlas/grid/GridBuilder.h"

using namespace std;
using namespace eckit;
using namespace atlas::grid;
namespace fs = boost::filesystem;


static void test_grids_from_grib_sample_directory( const std::string& directory);
static std::string determine_grib_samples_dir();

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

   grib_include_dir.replace(pos,grib_include_dir.length(),"/share/grib_api/samples");
   //std::cout << "GRIB_API_INCLUDE=" << grib_include_dir << "\n";

   return grib_include_dir;
}


static void test_grids_from_grib_sample_directory(const std::string& directory)
{
   fs::path full_path( fs::initial_path<fs::path>() );
   full_path = fs::system_complete( fs::path( directory ) );

   BOOST_CHECK(fs::exists( full_path ));
   BOOST_CHECK(fs::is_directory( full_path ));

   //std::cout << "\nIn directory: " << full_path.directory_string() << "\n\n";
   fs::directory_iterator end_iter;
   for ( fs::directory_iterator dir_itr( full_path ); dir_itr != end_iter; ++dir_itr )
   {
      try
      {
         fs::path relPath(directory + "/" + dir_itr->path().filename().string());

         // recurse down directories
         if ( fs::is_directory(dir_itr->status()) )  {
            test_grids_from_grib_sample_directory(relPath.string());
            continue;
         }

         std::cout << "\n==========================================================================" << std::endl;
         std::cout << "Opening GRIB file " << relPath.string() << std::endl;
         FILE* fp = fopen(relPath.string().c_str(),"r");
         BOOST_REQUIRE_MESSAGE(fp,"Could not open file " << relPath.string());

         std::cout << " Create a grib handle"<< std::endl;
         int err;
         grib_handle* handle = grib_handle_new_from_file(0,fp,&err);
         BOOST_REQUIRE_MESSAGE(err == 0,"grib_handle_new_from_file error " << err << " for file " << relPath.string());


         std::cout << " Get the grid type" << std::endl;
         char string_value[64];
         size_t len = sizeof(string_value)/sizeof(char);
         err = grib_get_string(handle,"gridType",string_value,&len);
         if ( err !=0 ) {
            std::cout << " grib_get_string(gridType) failed for file " << relPath.string() << " IGNORING !!!!\n";
            BOOST_REQUIRE_MESSAGE(fclose(fp) != -1,"error closing file " << relPath.string());
            continue;
         }
         std::string gridType = string_value;
         std::cout << " Create Grid derivatives " << gridType << std::endl;

         if ( gridType == "polar_stereographic" || gridType == "rotated_ll" || gridType == "reduced_ll" || gridType == "sh" ) {

            std::cout << " Ignoring grid types [ polar_stereographic | rotated_ll | reduced_ll || sh ] " << std::endl;
            std::cout << " close the grib file" << std::endl;
            err = grib_handle_delete(handle);
            BOOST_CHECK_MESSAGE(err == 0,"grib_handle_delete failed for " << relPath.string());

            // Close the file
            BOOST_REQUIRE_MESSAGE(fclose(fp) != -1,"error closing file " << relPath.string());
            continue;
         }


         // Unstructured grid can not handle Spherical harmonics
         atlas::grid::Grid::Ptr the_grid = GribGridBuilder::instance().build_grid_from_grib_handle(handle);
         BOOST_CHECK_MESSAGE(the_grid,"GribGridBuilder::instance().build_grid_from_grib_handle failed for file " << relPath.string());
         if (!the_grid) {
            BOOST_REQUIRE_MESSAGE(fclose(fp) != -1,"error closing file " << relPath.string());
            continue;
         }

         std::cout << " Check grid type are correct, found(" << gridType << ")" << std::endl;
         BOOST_CHECK_MESSAGE(the_grid->gridType() == gridType,"gridType(" << gridType << ") dir not match Grid constructor(" << the_grid->gridType() << ") for file " << relPath.string());

         std::cout << " close the grib file" << std::endl;
         err = grib_handle_delete(handle);
         BOOST_CHECK_MESSAGE(err == 0,"grib_handle_delete failed for " << relPath.string());

         // Close the file
         BOOST_REQUIRE_MESSAGE(fclose(fp) != -1,"error closing file " << relPath.string());
      }
      catch ( const std::exception & ex )
      {
         std::cout << dir_itr->path().filename() << " " << ex.what() << std::endl;
      }
   }
}
