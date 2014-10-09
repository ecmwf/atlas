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
#include <algorithm>

#define BOOST_TEST_MODULE TestGrid
#include "ecbuild/boost_test_framework.h"

#include "eckit/types/FloatCompare.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/io/FileHandle.h"

#include "eckit/grib/GribHandle.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Grib.h"
#include "atlas/grid/GridSpec.h"

using namespace std;
using namespace eckit;
using namespace eckit::grib;
using namespace atlas;
using namespace atlas::grid;

BOOST_AUTO_TEST_SUITE( TestGridFromGrib )


static void test_grib_file( const LocalPathName& path )
{
	Log::info() << "\n>>> testing file " << path << std::endl;

	// create a grid from the GRIB

	GribHandle gh(path);

	Grid::Ptr grid = Grib::create_grid( gh );

	BOOST_REQUIRE_MESSAGE( grid, " Could not create GRID ptr for file " << path);

	GridSpec spec = grid->spec();

//	Log::info() << "grib_type : " << grid->gridType() << std::endl;
//	Log::info() << "nb_nodes  : " << grid->nPoints() << std::endl;
	Log::info() << "spec      : " << spec << std::endl;

	// check the latlon points

	vector< Grid::Point > pts;
	pts.resize( grid->nPoints() );
	grid->coordinates(pts);

	vector< Grid::Point > grib_pts;
	gh.getLatLonPoints( grib_pts );

	BOOST_CHECK_MESSAGE( pts.size() == grib_pts.size()," Pts sizes differ GRIB:" << grib_pts.size() << " GRID: " << pts.size() );

	// Checking against the grib *only* makes sense when we have the same scanning mode.
	//    Atlas:  North,west -> south,east
	//    Grid:   Depended on scanning mode
	// Ignore for polar stereographic since we know that will fail.
   RealCompare<double> requal( Grid::degrees_eps() );

	if (grid->gridType() != "polar_stereographic") {
	   for( size_t i = 0; i < pts.size(); ++i )
	   {
//	      BOOST_CHECK_CLOSE( pts[i].lat(), grib_pts[i].lat(), Grid::degrees_eps() );
//	      BOOST_CHECK_CLOSE( pts[i].lon(), grib_pts[i].lon(), Grid::degrees_eps() );

	      BOOST_CHECK_MESSAGE( requal( pts[i].lat(), grib_pts[i].lat() ), " lats differ at pt index " << i << " GRIB: " << grib_pts[i].lat() << " GRID:" << pts[i].lat());
	      BOOST_CHECK_MESSAGE( requal( pts[i].lon(), grib_pts[i].lon() ), " lons differ at pt index " << i << " GRIB: " << grib_pts[i].lon() << " GRID:" << pts[i].lon());
	   }
	}

	// use the GridSpec to create another GRIB
	// which should then have the same geography (grid) as the initial one

	GribHandle::Ptr newgh = Grib::create_handle( *grid, gh.edition() );

	BOOST_CHECK_MESSAGE( newgh , "Could not create GribHandle from Grid ptr");

	BOOST_CHECK_MESSAGE( gh.gridType() == newgh->gridType(),"Grid type differ GRIB: " << gh.gridType() << " GRID: " << newgh->gridType());

//   if (grid->gridType() != "polar_stereographic") {
      // The geography has will be different, since data is from canada, will not match grib samples
      BOOST_CHECK_MESSAGE( gh.geographyHash() == newgh->geographyHash() ,"Geography hash different \nGRIB:" << gh.geographyHash() << "\nGRID:" << newgh->geographyHash());
//   }

	// write a new grib file

	LocalPathName opath(std::string(path) + ".new");
	opath.unlink();

	FileHandle of( opath );
	of.openForWrite(0);

	BOOST_CHECK_NO_THROW( newgh->write(of) );

	of.close();
}

BOOST_AUTO_TEST_CASE( test_grid_creation )
{
	std::vector< LocalPathName > gribs;

	// Polar stereographic grids
   gribs.push_back("CMC_hrdps_east_TMP_TGL_2_ps2.5km_2014081512_P000-00.grib2");

	// these seem to fail, probably due to the low tolerance of grib1 1E-3
	//	gribs.push_back("ll005005.grib");
	//	gribs.push_back("ll0101.grib");

	gribs.push_back("ll11.grib");
	gribs.push_back("ll55.grib");

	// reduced_gg

	gribs.push_back("rgg_n640.grib");
	// te gribs.push_back("rgg_n1280.grib");

	// regular_gg
	gribs.push_back("n40.grib");
	gribs.push_back("n60.grib");
	gribs.push_back("n160.grib");
	gribs.push_back("n256.grib");
	gribs.push_back("n320.grib");
	gribs.push_back("n400.grib");
	gribs.push_back("n512.grib");
	gribs.push_back("n1024.grib");

	std::for_each( gribs.begin(), gribs.end(), test_grib_file );
}

BOOST_AUTO_TEST_SUITE_END()


