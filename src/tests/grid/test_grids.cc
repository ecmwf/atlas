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

#include "atlas/library/Library.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Factory.h"
#include "eckit/types/FloatCompare.h"

#include "tests/AtlasTestEnvironment.h"

using StructuredGrid = atlas::grid::StructuredGrid;
using Grid       = atlas::Grid;
using Regular    = atlas::grid::RegularGrid;
using ReducedGaussianGrid    = atlas::grid::ReducedGaussianGrid;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_factory" )
{
  StructuredGrid structured = Grid("N80");

  Grid grid = Grid("N24");

  std::cout << "structured.ny() = " << structured.ny() << std::endl;
  std::cout << "grid.npts() = " << grid.size() << std::endl;

}

CASE( "test_regular_gg" )
{
  Regular grid( "F32" );

  EXPECT(grid.ny() == 64);
  EXPECT(grid.size() == 8192);
  // EXPECT(grid.type() == "regular_gaussian");

  // Full Gaussian Grid

  Grid::Config config;
  config.set("type","regular_gaussian");
  config.set("N",32);
  grid = Grid(config);
  EXPECT(grid.size() == 8192);
  // EXPECT(grid.type() == "regular_gaussian");



}

CASE( "test_reduced_gg" )
{
  StructuredGrid grid;

  grid = Grid( "N32" );
  EXPECT(grid.ny() == 64);
  EXPECT(grid.size() == 6114);

  grid = grid::ReducedGaussianGrid( {4,6,8,8,6,4} );

  EXPECT(grid.ny() == 6);
  EXPECT(grid.size() == 8+12+16);
}

CASE( "test_reduced_gg_ifs" )
{
  StructuredGrid grid( "N32" );

  // EXPECT(grid.N() ==    32);
  EXPECT(grid.ny() == 64);
  EXPECT(grid.size() == 6114);
  // EXPECT(grid.type() == "classic_gaussian");

}

CASE( "test_regular_ll" )
{
  // Constructor for N=8
  size_t nlon = 32;
  size_t nlat = 16;
  std::stringstream name; name << "Slat" << nlon << "x" << nlat;
  Regular grid( name.str() );

  EXPECT(grid.nx() == nlon);
  EXPECT(grid.ny() == nlat);
  EXPECT(grid.size() == 512);
  // EXPECT(grid.type() == "shifted_lat");
  EXPECT(grid.y(0) == 90.-0.5*(180./16.));
  EXPECT(grid.y(grid.ny()-1) == -90.+0.5*(180./16.));
  EXPECT(grid.x(0) == 0.);
  EXPECT(grid.x(grid.nx()-1) == 360.-360./32.);

  // Construct using builders/factories

  // Global Grid
  Grid::Config config1;
  config1.set("type","shifted_lat");
  config1.set("nx",32);
  config1.set("ny",16);
  grid = Grid(config1);
  EXPECT(grid.size() == 512);
  // EXPECT(gridptr->type() == "shifted_lat");

  Grid::Config config2;
  config2.set("type","shifted_lat");
  config2.set("N",8);
  grid = Grid(config2);
  EXPECT(grid.size() == 512);
  // EXPECT(gridptr->type() == "shifted_lat");

  Regular ll_poles( "L4x3" );
  EXPECT( ll_poles.nx() == 4);
  EXPECT( ll_poles.ny() == 3);

  Regular ll_nopoles( "Slat4x2" );
  EXPECT( ll_nopoles.nx() == 4);
  EXPECT( ll_nopoles.ny() == 2);
  EXPECT( eckit::types::is_approximately_equal( ll_nopoles.y(0), 45.) ); // tolerance was previously 1.e-5
  EXPECT( eckit::types::is_approximately_equal( ll_nopoles.y(1), -45. ) ); // tolerance was previously 1.e-5
  EXPECT( eckit::types::is_approximately_equal( ll_nopoles.x(0), 0. ) ); // tolerance was previously 1.e-5
  EXPECT( eckit::types::is_approximately_equal( ll_nopoles.x(1), 90. ) ); // tolerance was previously 1.e-5

}

CASE( "test_reducedgaussian" )
{
  StructuredGrid N640( "N640" );
  EXPECT(N640.size() == 2140702);
  ReducedGaussianGrid custom( N640.nx() );
  EXPECT(N640.size() == custom.size());
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    return atlas::test::run( argc, argv );
}
