/*
 * (C) Copyright 1996-2016 ECMWF.
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

#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"

#include "atlas/atlas.h"
#include "atlas/grid/Grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/grid/grids.h"


using namespace eckit;
using namespace atlas;
using namespace atlas::grid;

namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( init ) { eckit::mpi::init(); atlas::grid::load(); }

BOOST_AUTO_TEST_CASE( test_factory )
{
  SharedPtr<global::Structured> reducedgrid( global::Structured::create("N80") );

  SharedPtr<Grid> grid ( Grid::create("N24") );

  std::cout << "reducedgrid->nlat() = " << reducedgrid->nlat() << std::endl;
  std::cout << "grid->npts() = " << grid->npts() << std::endl;

}

BOOST_AUTO_TEST_CASE( test_regular_gg )
{
  // Constructor for N=32
  global::gaussian::RegularGaussian grid(32);

  BOOST_CHECK_EQUAL(grid.N(), 32);
  BOOST_CHECK_EQUAL(grid.nlat(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 8192);
  BOOST_CHECK_EQUAL(grid.gridType(),"regular_gaussian");

  // Construct using builders/factories

  Grid::Ptr gridptr;

  // Full Gaussian Grid

  util::Config spec;
  spec.set("grid_type","regular_gaussian");
  spec.set("N",32);
  gridptr = Grid::Ptr( Grid::create(spec) );
  BOOST_CHECK_EQUAL(gridptr->npts(), 8192);
  BOOST_CHECK_EQUAL(gridptr->gridType(),"regular_gaussian");

  // Add bounding box to spec
  BoundBox bbox( 90., 0., 180., 0.);
  spec.set("bbox_s", bbox.min().lat());
  spec.set("bbox_w", bbox.min().lon());
  spec.set("bbox_n", bbox.max().lat());
  spec.set("bbox_e", bbox.max().lon());
  gridptr = Grid::Ptr( Grid::create(spec) );

}


BOOST_AUTO_TEST_CASE( test_reduced_gg )
{
  long nlon[] = {4,6,8};
  grid::global::gaussian::ReducedGaussian grid(3,nlon);
  BOOST_CHECK_EQUAL(grid.N(),3);
  BOOST_CHECK_EQUAL(grid.nlat(),6);
  BOOST_CHECK_EQUAL(grid.npts(),8+12+16);
  BOOST_CHECK_EQUAL(grid.gridType(),"reduced_gaussian");
}

BOOST_AUTO_TEST_CASE( test_reduced_gg_ifs )
{
  grid::global::gaussian::ClassicGaussian grid(32);

  BOOST_CHECK_EQUAL(grid.N(),    32);
  BOOST_CHECK_EQUAL(grid.nlat(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 6114);
  BOOST_CHECK_EQUAL(grid.gridType(),"classic_gaussian");

}

BOOST_AUTO_TEST_CASE( test_regular_ll )
{
  // Constructor for N=8
  size_t nlon = 32;
  size_t nlat = 16;
  global::lonlat::ShiftedLat grid(nlon,nlat);

  BOOST_CHECK_EQUAL(grid.nlon(), nlon);
  BOOST_CHECK_EQUAL(grid.nlat(), nlat);
  BOOST_CHECK_EQUAL(grid.npts(), 512);
  BOOST_CHECK_EQUAL(grid.gridType(),"shifted_lat");
  BOOST_CHECK_EQUAL(grid.lat(0), 90.-0.5*(180./16.));
  BOOST_CHECK_EQUAL(grid.lat(grid.nlat()-1), -90.+0.5*(180./16.));
  BOOST_CHECK_EQUAL(grid.lon(0), 0.);
  BOOST_CHECK_EQUAL(grid.lon(grid.nlon()-1), 360.-360./32.);


  // Construct using builders/factories

  Grid::Ptr gridptr;

  // Global Grid
  util::Config spec;
  spec.set("grid_type","shifted_lat");
  spec.set("nlon",32);
  spec.set("nlat",16);
  gridptr = Grid::Ptr( Grid::create(spec) );
  BOOST_CHECK_EQUAL(gridptr->npts(), 512);
  BOOST_CHECK_EQUAL(gridptr->gridType(),"shifted_lat");

  util::Config spec2;
  spec2.set("grid_type","shifted_lat");
  spec2.set("N",8);
  gridptr = Grid::Ptr( Grid::create(spec2) );
  BOOST_CHECK_EQUAL(gridptr->npts(), 512);
  BOOST_CHECK_EQUAL(gridptr->gridType(),"shifted_lat");

  global::lonlat::RegularLonLat ll_poles(90.,90.);
  BOOST_CHECK_EQUAL( ll_poles.nlat(), 3);
  BOOST_CHECK_EQUAL( ll_poles.nlon(), 4);

  global::lonlat::ShiftedLat ll_nopoles(90.,90.);
  BOOST_CHECK_EQUAL( ll_nopoles.nlat(), 2);
  BOOST_CHECK_EQUAL( ll_nopoles.nlon(), 4);
  BOOST_CHECK_CLOSE( ll_nopoles.lat(0), 45., 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.lat(1), -45. , 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.lon(0), 0. , 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.lon(1), 90. , 1.e-5);
}


BOOST_AUTO_TEST_CASE( finalize ) { eckit::mpi::finalize(); }

} // namespace test
} // namespace atlas
