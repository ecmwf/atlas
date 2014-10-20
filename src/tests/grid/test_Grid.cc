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
#include "ecbuild/boost_test_framework.h"

#include "eckit/io/StdFile.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/grib/GribHandle.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Grib.h"
#include "atlas/grid/GridSpec.h"

using namespace std;
using namespace eckit;
using namespace eckit::grib;
using namespace atlas;
using namespace atlas::grid;

/// Test for Grid* derivatives
/// This test uses the grib samples directory.
/// We open each file and attempt to create the Grid* class.
///
/// For the full test run with ctest -VV

static void test_grib_file(const std::string& file);

BOOST_AUTO_TEST_SUITE( TestGrid )

struct ArgsFixture {
   ArgsFixture(): argc(boost::framework::master_test_suite().argc),
           argv(boost::framework::master_test_suite().argv){}
   int argc;
   char **argv;
};

BOOST_FIXTURE_TEST_CASE ( test_Grid, ArgsFixture ) {
    cout << "Grid:: ...test_Grid argc = " << argc << " ";
    if (argc == 2) cout << argv[1];
    cout << "\n";
    BOOST_REQUIRE_MESSAGE( argc == 2, "You miss one argument" );
    BOOST_REQUIRE_MESSAGE( argv[1] != "some_required_arg", "The first arg it's wrong!!");

    test_grib_file(argv[1]);
}

BOOST_AUTO_TEST_SUITE_END()


static void test_grib_file(const std::string& fpath)
{
	std::cout << std::endl;

	cout << "-----------------------------------------------------------------------" << endl;
	std::cout << "GRIB file " << fpath << std::endl;

	LocalPathName path(fpath);

	eckit::grib::GribHandle gh( path );

	std::string gridType = gh.gridType();

	std::cout << "gridType : " << gridType << std::endl;

	// skip polar_stereographic
	if ( gridType == "polar_stereographic" || gridType == "sh" || gridType.empty() )
	{
		std::cout << " ** Ignoring grid types [ polar_stereographic | sh | "" ] " << std::endl;
		return;
	}

	long editionNumber = gh.edition();

	atlas::grid::Grid::Ptr g = Grib::create_grid(gh);

	BOOST_CHECK_MESSAGE(g,"Grib::create_grid failed for file " << fpath);
	if (!g) return;

	// The Grid produced, has a GRID spec, the grid spec can be used to,
	// make sure the grid types match
	GridSpec g_spec = g->spec();
	BOOST_CHECK_MESSAGE(g->gridType() == gridType,"gridType(" << gridType << ") did not match Grid constructor(" << g->gridType() << ") for file " << fpath);
	BOOST_CHECK_MESSAGE(g_spec.grid_type() == gridType,"gridType(" << gridType << ") did not match GridSpec constructor(" << g_spec.grid_type() << ") for file " << fpath);

//	// find the corresponding sample file.
//	// However we need to take into account that the GRIB samples, file are *NOT* unique in their GRID definition.
//	// The sample file name produced does not have '.tmpl' extension
//	std::string generated_sample_file_name = Grib::grib_sample_file( g_spec , editionNumber);
//	BOOST_CHECK_MESSAGE( !generated_sample_file_name.empty()," Could *not* find sample file for grid_spec " << g_spec );
//
//
//	// Note: many of the grib samples files are not UNIQUE in their grid specification:
//	// hence the use of WARN.
//	// remove .tmpl and get base part
//
//	LocalPathName base_name = path.baseName(false);
//	std::string grib_sample_file = base_name.localPath();
//	BOOST_WARN_MESSAGE( generated_sample_file_name == grib_sample_file, "\nCould not match samples expected '"
//						<< grib_sample_file << "' but found('"
//						<< generated_sample_file_name
//						<< "') for grid spec "
//						<< g_spec );
}
