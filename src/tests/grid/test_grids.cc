/*
 * (C) Copyright 1996-2017 ECMWF.
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


#include "atlas/grid/detail/grid/regular/RegularGaussian.h"
#include "atlas/grid/detail/grid/regular/RegularLonLat.h"
#include "atlas/grid/detail/grid/regular/ShiftedLat.h"
#include "atlas/grid/detail/grid/reduced/ClassicGaussian.h"


#include "atlas/atlas.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Factory.h"

using Structured = atlas::grid::StructuredGrid;
using Grid       = atlas::grid::Grid;
using Regular    = atlas::grid::RegularGrid;


namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( test_factory )
{
  Structured structured = Grid("N80");

  Grid grid = Grid("N24");

  std::cout << "structured.ny() = " << structured.ny() << std::endl;
  std::cout << "grid.npts() = " << grid.npts() << std::endl;

}

BOOST_AUTO_TEST_CASE( test_regular_gg )
{
  // Constructor for N=32
  Regular grid( new grid::detail::grid::regular::RegularGaussian(32) );

  BOOST_CHECK_EQUAL(grid.ny(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 8192);
  // BOOST_CHECK_EQUAL(grid.gridType(),"regular_gaussian");

  // Full Gaussian Grid

  Grid::Config config;
  config.set("type","regular_gaussian");
  config.set("N",32);
  grid = Grid(config);
  BOOST_CHECK_EQUAL(grid.npts(), 8192);
  // BOOST_CHECK_EQUAL(grid.gridType(),"regular_gaussian");



}

BOOST_AUTO_TEST_CASE( test_reduced_gg )
{
  Structured grid;
  
  grid = Grid( new grid::detail::grid::reduced::ClassicGaussian(32) );
  BOOST_CHECK_EQUAL(grid.ny(),64);
  BOOST_CHECK_EQUAL(grid.npts(),6114);

  long nlon[] = {4,6,8};
  grid = Grid( new grid::detail::grid::reduced::ReducedGaussian(3,nlon) );
  BOOST_CHECK_EQUAL(grid.ny(),6);
  BOOST_CHECK_EQUAL(grid.npts(),8+12+16);
}

BOOST_AUTO_TEST_CASE( test_reduced_gg_ifs )
{
  Structured grid( new grid::detail::grid::reduced::ClassicGaussian(32) );

  // BOOST_CHECK_EQUAL(grid.N(),    32);
  BOOST_CHECK_EQUAL(grid.ny(), 64);
  BOOST_CHECK_EQUAL(grid.npts(), 6114);
  // BOOST_CHECK_EQUAL(grid.gridType(),"classic_gaussian");

}

BOOST_AUTO_TEST_CASE( test_regular_ll )
{
  // Constructor for N=8
  size_t nlon = 32;
  size_t nlat = 16;
  Regular grid( new grid::detail::grid::regular::ShiftedLat(nlon,nlat) );

  BOOST_CHECK_EQUAL(grid.nx(), nlon);
  BOOST_CHECK_EQUAL(grid.ny(), nlat);
  BOOST_CHECK_EQUAL(grid.npts(), 512);
  // BOOST_CHECK_EQUAL(grid.gridType(),"shifted_lat");
  BOOST_CHECK_EQUAL(grid.y(0), 90.-0.5*(180./16.));
  BOOST_CHECK_EQUAL(grid.y(grid.ny()-1), -90.+0.5*(180./16.));
  BOOST_CHECK_EQUAL(grid.x(0), 0.);
  BOOST_CHECK_EQUAL(grid.x(grid.nx()-1), 360.-360./32.);

  // Construct using builders/factories

  // Global Grid
  Grid::Config config1;
  config1.set("type","shifted_lat");
  config1.set("nx",32);
  config1.set("ny",16);
  grid = Grid(config1);
  BOOST_CHECK_EQUAL(grid.npts(), 512);
  // BOOST_CHECK_EQUAL(gridptr->gridType(),"shifted_lat");

  Grid::Config config2;
  config2.set("type","shifted_lat");
  config2.set("N",8);
  grid = Grid(config2);
  BOOST_CHECK_EQUAL(grid.npts(), 512);
  // BOOST_CHECK_EQUAL(gridptr->gridType(),"shifted_lat");

  Regular ll_poles( new grid::detail::grid::regular::RegularLonLat(4, 3) );
  BOOST_CHECK_EQUAL( ll_poles.nx(), 4);
  BOOST_CHECK_EQUAL( ll_poles.ny(), 3);

  Regular ll_nopoles( new grid::detail::grid::regular::ShiftedLat(4, 2) );
  BOOST_CHECK_EQUAL( ll_nopoles.nx(), 4);
  BOOST_CHECK_EQUAL( ll_nopoles.ny(), 2);
  BOOST_CHECK_CLOSE( ll_nopoles.y(0), 45., 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.y(1), -45. , 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.x(0), 0. , 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.x(1), 90. , 1.e-5);
  
}

BOOST_AUTO_TEST_CASE( test_reducedgaussian )
{
  Structured N640( new grid::detail::grid::reduced::ClassicGaussian(640) );
  BOOST_CHECK_EQUAL(N640.npts(),2140702);
  Structured custom( new grid::detail::grid::reduced::ReducedGaussian(N640.ny()/2, N640.nx().data()) );
  BOOST_CHECK_EQUAL(N640.npts(),custom.npts());
}

} // namespace test
} // namespace atlas
