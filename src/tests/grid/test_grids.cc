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

#include "atlas/library/Library.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Factory.h"

using StructuredGrid = atlas::grid::StructuredGrid;
using Grid       = atlas::Grid;
using Regular    = atlas::grid::RegularGrid;
using ReducedGaussianGrid    = atlas::grid::ReducedGaussianGrid;


namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( test_factory )
{
  StructuredGrid structured = Grid("N80");

  Grid grid = Grid("N24");

  std::cout << "structured.ny() = " << structured.ny() << std::endl;
  std::cout << "grid.npts() = " << grid.size() << std::endl;

}

BOOST_AUTO_TEST_CASE( test_regular_gg )
{
  Regular grid( "F32" );

  BOOST_CHECK_EQUAL(grid.ny(), 64);
  BOOST_CHECK_EQUAL(grid.size(), 8192);
  // BOOST_CHECK_EQUAL(grid.type(),"regular_gaussian");

  // Full Gaussian Grid

  Grid::Config config;
  config.set("type","regular_gaussian");
  config.set("N",32);
  grid = Grid(config);
  BOOST_CHECK_EQUAL(grid.size(), 8192);
  // BOOST_CHECK_EQUAL(grid.type(),"regular_gaussian");



}

BOOST_AUTO_TEST_CASE( test_reduced_gg )
{
  StructuredGrid grid;

  grid = Grid( "N32" );
  BOOST_CHECK_EQUAL(grid.ny(),64);
  BOOST_CHECK_EQUAL(grid.size(),6114);

  grid = grid::ReducedGaussianGrid( {4,6,8,8,6,4} );

  BOOST_CHECK_EQUAL(grid.ny(),6);
  BOOST_CHECK_EQUAL(grid.size(),8+12+16);
}

BOOST_AUTO_TEST_CASE( test_reduced_gg_ifs )
{
  StructuredGrid grid( "N32" );

  // BOOST_CHECK_EQUAL(grid.N(),    32);
  BOOST_CHECK_EQUAL(grid.ny(), 64);
  BOOST_CHECK_EQUAL(grid.size(), 6114);
  // BOOST_CHECK_EQUAL(grid.type(),"classic_gaussian");

}

BOOST_AUTO_TEST_CASE( test_regular_ll )
{
  // Constructor for N=8
  size_t nlon = 32;
  size_t nlat = 16;
  std::stringstream name; name << "Slat" << nlon << "x" << nlat;
  Regular grid( name.str() );

  BOOST_CHECK_EQUAL(grid.nx(), nlon);
  BOOST_CHECK_EQUAL(grid.ny(), nlat);
  BOOST_CHECK_EQUAL(grid.size(), 512);
  // BOOST_CHECK_EQUAL(grid.type(),"shifted_lat");
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
  BOOST_CHECK_EQUAL(grid.size(), 512);
  // BOOST_CHECK_EQUAL(gridptr->type(),"shifted_lat");

  Grid::Config config2;
  config2.set("type","shifted_lat");
  config2.set("N",8);
  grid = Grid(config2);
  BOOST_CHECK_EQUAL(grid.size(), 512);
  // BOOST_CHECK_EQUAL(gridptr->type(),"shifted_lat");

  Regular ll_poles( "L4x3" );
  BOOST_CHECK_EQUAL( ll_poles.nx(), 4);
  BOOST_CHECK_EQUAL( ll_poles.ny(), 3);

  Regular ll_nopoles( "Slat4x2" );
  BOOST_CHECK_EQUAL( ll_nopoles.nx(), 4);
  BOOST_CHECK_EQUAL( ll_nopoles.ny(), 2);
  BOOST_CHECK_CLOSE( ll_nopoles.y(0), 45., 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.y(1), -45. , 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.x(0), 0. , 1.e-5);
  BOOST_CHECK_CLOSE( ll_nopoles.x(1), 90. , 1.e-5);

}

BOOST_AUTO_TEST_CASE( test_reducedgaussian )
{
  StructuredGrid N640( "N640" );
  BOOST_CHECK_EQUAL(N640.size(),2140702);
  ReducedGaussianGrid custom( N640.nx() );
  BOOST_CHECK_EQUAL(N640.size(),custom.size());
}

} // namespace test
} // namespace atlas
