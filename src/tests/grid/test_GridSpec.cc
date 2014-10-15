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

#define BOOST_TEST_MODULE TestGridSpec
#include "ecbuild/boost_test_framework.h"

#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"

#include "eckit/grib/GribHandle.h"
#include "eckit/grib/GribAccessor.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Grib.h"
#include "atlas/grid/GridSpec.h"
#include "atlas/grid/ReducedGG.h"


using namespace std;
using namespace eckit;
using namespace eckit::grib;
using namespace atlas;
using namespace atlas::grid;

/// Test for Grid* derivatives
/// This test uses the grib samples directory.
/// We open each file and attempt to create the Grid* class.
/// However this has exposed errors in the grib samples files,
/// especially for reduced gaussian grids. Hence we will need to wait till Shahram has
/// has a chance to fix the issues, before all the tests pass.

static void test_grids_from_grib_sample_directory( const std::string& directory);
static void test_grib_file(const std::string& file);
static void read_data_points(grib_handle *h, std::vector<Grid::Point> & points);
static void comparePointList(const std::vector<Grid::Point>& grib_pntlist, const std::vector<Grid::Point>& points, double epsilon, eckit::grib::GribHandle& gh);


BOOST_AUTO_TEST_SUITE( TestGridSpec )

BOOST_AUTO_TEST_CASE( test_grib_to_grid_to_gridspec )
{
   cout << "Grid:: ...test_grib_to_grid_to_gridspec\n";

   // Traverse all the GRIB samples files, for gridType first determine sample dir
   std::vector<std::string> sample_dirs;
   Grib::determine_grib_samples_dir(sample_dirs);
   BOOST_REQUIRE_MESSAGE(!sample_dirs.empty(),"Expected sample dirs to be found");

   // now test these dirs
   for(size_t i = 0; i < sample_dirs.size(); ++i) {
      test_grids_from_grib_sample_directory(sample_dirs[i]);
   }
}


BOOST_AUTO_TEST_CASE( test_rotated_grids )
{
   cout << "Grid:: ...test_rotated_grids \n";

   // Note: we need to wait till grib iterator, rotates the points.
   // At the moment(grib 13.1) it just, return the points, in regular lat long fashion.
   std::string path = "/scratch/ma/ma0/wind_rotated_latlon.grb";
   test_grib_file( path );
}

BOOST_AUTO_TEST_SUITE_END()


static void test_grids_from_grib_sample_directory(const std::string& directory)
{
   cout << "*********************************************************************************\n";
   cout << "traversing directory " << directory << "\n";
   PathName dir_path(directory);
   BOOST_CHECK(dir_path.exists());
   BOOST_CHECK(dir_path.isDir());

   //int count = 0;
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

static void test_grib_file(const std::string& fpath)
{
   std::cout << "\n===============================================================================================================" << std::endl;
   LocalPathName path(fpath);
   std::cout << "Opening GRIB file " << fpath << std::endl;
   LocalPathName baseName = path.baseName(false);
   if (baseName.path() == "budg") {
      std::cout << "Ignoring " << path.baseName() << " not a grid based grib file\n";
      return;
   }

   eckit::grib::GribHandle gh(path);
   std::string gridType = gh.gridType();
   std::cout << " Create Grid derivatives " << gridType << std::endl;

   if ( gridType == "polar_stereographic" || gridType == "sh" || gridType.empty())
   {
      // The GRIB samples for polar stereographic are corrupt in terms of the GRID section, hence IGNORE
      std::cout << " ** Ignoring grid types [ polar_stereographic | sh ] " << std::endl;
      return;
   }

   // Create Grid derivatives from the GRIB file
   atlas::grid::Grid::Ptr grid_created_from_grib = Grib::create_grid(gh);
   BOOST_CHECK_MESSAGE(grid_created_from_grib,"GRIBGridBuilder::instance().build_grid_from_grib_handle failed for file " << fpath);
   if (!grid_created_from_grib) return;


   // The Grid produced, has a GRID spec, the grid spec can be used to, make sure the grid types match
   GridSpec g_spec = grid_created_from_grib->spec();
   std::cout << " " << g_spec << std::endl;
   BOOST_CHECK_MESSAGE(grid_created_from_grib->gridType() == gridType,"gridType(" << gridType << ") did not match Grid constructor(" << grid_created_from_grib->gridType() << ") for file " << fpath);
   BOOST_CHECK_MESSAGE(g_spec.grid_type() == gridType,"gridType(" << gridType << ") did not match GridSpec constructor(" << g_spec.grid_type() << ") for file " << fpath);


   // From the Spec, create another Grid, we should get back the same Grid
   Grid::Ptr grid_created_from_spec = Grid::create(g_spec);
   BOOST_CHECK_MESSAGE(grid_created_from_spec,"Failed to create GRID from GridSpec");
   bool grid_compare = grid_created_from_grib->same(*grid_created_from_spec);
   BOOST_CHECK_MESSAGE(grid_compare,"The grids are different");


   // For reduced Guassian Grid, check the no of pts per latitude read from grib, matches the computed values
   if (gridType == "reduced_gg") {

      ReducedGG* read_from_grib = dynamic_cast<ReducedGG*>(grid_created_from_grib.get());
      BOOST_CHECK_MESSAGE(read_from_grib,"Downcast to ReducedGG failed ?");

      if (read_from_grib) {
         // Create on the fly, this will compute no of pts per latitude on the fly
         ReducedGG reducedgg(read_from_grib->gaussianNumber());
         BOOST_CHECK_MESSAGE(read_from_grib->pointsPerLatitude() == reducedgg.pointsPerLatitude(),"Pts per latitide read from grib, different to computed pts per latitude");
      }
   }


   // epsilon varies depending on the edition number
   long editionNumber = GribAccessor<long>("editionNumber")(gh);
   BOOST_CHECK_MESSAGE(editionNumber ==  gh.edition(),"Edition numbers dont match");

   double epsilon = (editionNumber == 1) ? 1e-3 : 1e-6;

   if (gridType == "reduced_gg" || gridType == "regular_gg") {

      // ---------------------------------------------------------------------------------------
      // Grib samples file have errors in the EXPECTED_longitudeOfLastGridPointInDegrees
      // This can then effect, comparison of the points
      // ----------------------------------------------------------------------------------------
      double guass = GribAccessor<long>("numberOfParallelsBetweenAPoleAndTheEquator")(gh);
      Grid::BoundBox bbox = grid_created_from_grib->boundingBox();

      double EXPECTED_longitudeOfLastGridPointInDegrees = 360.0 - (90.0/guass);
      std::cout << " EXPECTED longitudeOfLastGridPointInDegrees     " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << EXPECTED_longitudeOfLastGridPointInDegrees << std::endl;
      std::cout << " east                                           " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << bbox.east() << std::endl;
      BOOST_CHECK_CLOSE(bbox.east(),EXPECTED_longitudeOfLastGridPointInDegrees,epsilon);
   }


   // get GRIB points
   std::vector<Grid::Point> grib_pntlist;
   gh.getLatLonPoints( grib_pntlist );
   BOOST_CHECK_MESSAGE( grid_created_from_grib->nPoints() == grib_pntlist.size(),"GRIB pt list size " << grib_pntlist.size() << " different to GRID " << grid_created_from_grib->nPoints());

   // get the GRID points
   std::vector<Grid::Point> grid_points; grid_points.resize( grid_created_from_grib->nPoints());
   grid_created_from_grib->coordinates(grid_points);

   comparePointList(grib_pntlist,grid_points,epsilon,gh);
}

void comparePointList(const std::vector<Grid::Point>& grib_pntlist, const std::vector<Grid::Point>& points, double epsilon, eckit::grib::GribHandle& gh)
{
   BOOST_CHECK_MESSAGE(  points.size() == grib_pntlist.size(),"\n  **GRIB pt list size " << grib_pntlist.size() << " different to GRID " << points.size() );

   RealCompare<double> isEqual(epsilon);

   int print_point_list = 0;
   std::streamsize old_precision = cout.precision();
   for(size_t i =0;  i < points.size() && i < grib_pntlist.size(); i++) {
      if (!isEqual(points[i].lat(),grib_pntlist[i].lat()) ||
          !isEqual(points[i].lon(),grib_pntlist[i].lon()))
      {
         if (print_point_list == 0) {
            Log::info() << " Point list DIFFER, show first 10, epsilon(" << epsilon << ")" << std::endl;
            Log::info() << setw(40) << left << "     GRID   " << setw(40) << left << "         GRIB"<< std::endl;
         }
         Log::info() << setw(3) << i << " :"
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << points[i].lat() << ", "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << points[i].lon() << "  "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << grib_pntlist[i].lat() << ", "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << grib_pntlist[i].lon()
                  << std::endl;
         print_point_list++;
         if (print_point_list > 10) break;
      }
   }

   // reset precision
   Log::info() << std::setprecision(old_precision);
}

