/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <algorithm>
#include <iomanip>
#define BOOST_TEST_MODULE TestGrids
#include "ecbuild/boost_test_framework.h"

#include <eckit/memory/Factory.h>
#include <eckit/memory/Builder.h>
#include "atlas/atlas.h"
#include "atlas/Grid.h"
#include "atlas/GridSpec.h"
#include "atlas/GridSpecParams.h"
#include "atlas/grids/grids.h"

using namespace eckit;
using namespace atlas;
using namespace atlas::grids;

BOOST_AUTO_TEST_CASE( init ) { MPL::init(); }

BOOST_AUTO_TEST_CASE( test_factory )
{
  ReducedGrid::Ptr reducedgrid( ReducedGrid::create("reduced_gg.N80") );

  Grid::Ptr grid ( Grid::create("reduced_gg.N24") );

  std::cout << "reducedgrid->nlat() = " << reducedgrid->nlat() << std::endl;
  std::cout << "grid->npts() = " << grid->npts() << std::endl;

}

BOOST_AUTO_TEST_CASE( test_regular_gg )
{
  // Constructor for N=32
  grids::GaussianGrid grid(32);
  GridSpec spec = grid.spec();

  BOOST_CHECK_EQUAL(grid.N(), 32);
  BOOST_CHECK_EQUAL(grid.nlat(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 8192);
  BOOST_CHECK_EQUAL(grid.grid_type(),"regular_gg");

  // Cropping with global boundbox should not do anything
  grid.crop( Grid::BoundBox( 90., -90., 360.-grid.degrees_eps(), 0.).global(false) );

  BOOST_CHECK_EQUAL(grid.nlat(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 8192);

  // Crop
  grid.crop( Grid::BoundBox( 90., 0., 180., 0.).global(false) );

  BOOST_CHECK_EQUAL(grid.nlat(), 32);
  BOOST_CHECK_EQUAL(grid.npts(), 2080);

  // Construct using builders/factories

  Grid::Ptr gridptr;

  // Full Gaussian Grid
  gridptr = Grid::Ptr( Grid::create(GridSpecParams(spec)) );
  BOOST_CHECK_EQUAL(gridptr->npts(), 8192);
  BOOST_CHECK_EQUAL(gridptr->grid_type(),"regular_gg");
  gridptr = Grid::create(spec);
  BOOST_CHECK_EQUAL(gridptr->npts(), 8192);
  BOOST_CHECK_EQUAL(gridptr->grid_type(),"regular_gg");

  // Add bounding box to spec
  spec.set_bounding_box( Grid::BoundBox( 90., 0., 180., 0.) );
  gridptr = Grid::Ptr( Grid::create(spec) );
  BOOST_CHECK_EQUAL(gridptr->npts(), 2080);
  BOOST_CHECK_EQUAL(gridptr->grid_type(),"regular_gg");
  gridptr = Grid::Ptr( Grid::create(spec) );
  BOOST_CHECK_EQUAL(gridptr->npts(), 2080);
  BOOST_CHECK_EQUAL(gridptr->grid_type(),"regular_gg");

  GridSpec spec2("regular_gg");
  spec2.set("N",16);
  Grid::Ptr N16 ( Grid::create(spec2) );
  BOOST_CHECK_EQUAL(N16->npts(), 2048);
  BOOST_CHECK_EQUAL(N16->grid_type(),"regular_gg");
}


BOOST_AUTO_TEST_CASE( test_reduced_gg )
{
  int nlon[] = {4,6,8};
  grids::ReducedGaussianGrid grid(3,nlon);
  BOOST_CHECK_EQUAL(grid.N(),3);
  BOOST_CHECK_EQUAL(grid.nlat(),6);
  BOOST_CHECK_EQUAL(grid.npts(),8+12+16);
  BOOST_CHECK_EQUAL(grid.grid_type(),"reduced_gg");
}

BOOST_AUTO_TEST_CASE( test_reduced_gg_ifs )
{
  grids::reduced_gg::N32 grid;

  BOOST_CHECK_EQUAL(grid.N(),    32);
  BOOST_CHECK_EQUAL(grid.nlat(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 6114);
  BOOST_CHECK_EQUAL(grid.grid_type(),"reduced_gg");


  grid.crop( Grid::BoundBox( 90., -90., 360.-grid.degrees_eps(), 0.).global(false) );

  BOOST_CHECK_EQUAL(grid.nlat(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 6114);

  grid.crop( Grid::BoundBox( 90., 0., 180., 0.).global(false) );

  BOOST_CHECK_EQUAL(grid.nlat(), 32);
  BOOST_CHECK_EQUAL(grid.npts(), 1559);
}

BOOST_AUTO_TEST_CASE( finalize ) { MPL::finalize(); }
