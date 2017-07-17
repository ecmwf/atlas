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

#include "tests/AtlasFixture.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"


#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"

using Grid           = atlas::Grid;
using StructuredGrid = atlas::grid::StructuredGrid;
using RegularGrid    = atlas::grid::RegularGrid;
using Config         = atlas::util::Config;

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_from_string_L32 )
{
  Grid grid;
  BOOST_CHECK( not grid );

  grid = Grid("L32");
  BOOST_CHECK( grid );
  BOOST_CHECK( StructuredGrid(grid) == true  );
  BOOST_CHECK( RegularGrid(grid)    == true  );

  auto structured = StructuredGrid(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 65 );
  BOOST_CHECK_EQUAL( structured.periodic(), true );
  BOOST_CHECK_EQUAL( structured.nx(0), 128 );
  BOOST_CHECK_EQUAL( structured.y().front(), 90. );
  BOOST_CHECK_EQUAL( structured.y().back(), -90. );

  auto regular = RegularGrid(grid);
  BOOST_CHECK_EQUAL( regular.ny(), 65 );
  BOOST_CHECK_EQUAL( regular.periodic(), true );
  BOOST_CHECK_EQUAL( regular.nx(), 128 );
  BOOST_CHECK_EQUAL( regular.y().front(), 90. );
  BOOST_CHECK_EQUAL( regular.y().back(), -90. );

}

BOOST_AUTO_TEST_CASE( test_from_string_O32 )
{
  Grid grid;
  BOOST_CHECK( not grid );

  grid = Grid("O32");
  BOOST_CHECK( grid );

  BOOST_CHECK_EQUAL( StructuredGrid(grid), true  );
  BOOST_CHECK_EQUAL( RegularGrid(grid)   , false );

  auto structured = StructuredGrid(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 64 );
  BOOST_CHECK_EQUAL( structured.periodic(), true );
  BOOST_CHECK_EQUAL( structured.nx().front(), 20 );
}

BOOST_AUTO_TEST_CASE( test_from_string_O32_with_domain )
{
  Grid grid;
  BOOST_CHECK( not grid );

  grid = Grid("O32",RectangularDomain( {0,90}, {0,90} ) );
  BOOST_CHECK( grid );

  BOOST_CHECK_EQUAL( StructuredGrid(grid), true  );
  BOOST_CHECK_EQUAL( RegularGrid(grid)   , false );

  auto structured = StructuredGrid(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 32 );
  BOOST_CHECK_EQUAL( structured.periodic(), false );
  BOOST_CHECK_EQUAL( structured.nx().front(), 6 );

  output::Gmsh gmsh("test_grid_ptr_O32_subdomain.msh");
  Mesh mesh = meshgenerator::StructuredMeshGenerator().generate(grid);
  gmsh.write(mesh);

}


BOOST_AUTO_TEST_CASE( test_structured_1 )
{
  std::stringstream json;
  json <<
      "{"
          "\"type\" : \"structured\","
          "\"yspace\" : { \"type\":\"linear\", \"N\":9,  \"start\":90, \"end\":-90 },"
          "\"xspace\" : { \"type\":\"linear\", \"N\":16, \"start\":0,  \"end\":360, \"endpoint\":false }"
      "}";
  json.seekp(0);

  Grid grid;
  BOOST_CHECK( not grid );

  Config json_config(json);

  grid = StructuredGrid( json_config );
  BOOST_CHECK( grid );
  BOOST_CHECK( StructuredGrid(grid) == true  );
  BOOST_CHECK( RegularGrid(grid)    == true  );

  auto structured = StructuredGrid(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 9 );
  BOOST_CHECK_EQUAL( structured.periodic(), true );
  BOOST_CHECK_EQUAL( structured.nx(0), 16 );
  BOOST_CHECK_EQUAL( structured.y().front(), 90. );
  BOOST_CHECK_EQUAL( structured.y().back(), -90. );

  auto regular = RegularGrid(grid);
  BOOST_CHECK_EQUAL( regular.ny(), 9 );
  BOOST_CHECK_EQUAL( regular.periodic(), true );
  BOOST_CHECK_EQUAL( regular.nx(), 16 );
  BOOST_CHECK_EQUAL( regular.y().front(), 90. );
  BOOST_CHECK_EQUAL( regular.y().back(), -90. );

  output::Gmsh gmsh("test_grid_ptr.msh");
  Mesh mesh = meshgenerator::StructuredMeshGenerator().generate(grid);
  gmsh.write(mesh);
}


BOOST_AUTO_TEST_CASE( test_structured_2 )
{
  using XSpace = StructuredGrid::XSpace;
  using YSpace = StructuredGrid::YSpace;
  using Domain = StructuredGrid::Domain;
  using Projection = StructuredGrid::Projection;
  StructuredGrid grid(
      XSpace( {0.,360.} , {2,4,6,6,4,2} , false ),
      YSpace( grid::LinearSpacing( {90.,-90.}, 6 ) ),
      Projection(),
      Domain() );
  BOOST_CHECK( grid );

  output::Gmsh gmsh("test_grid_ptr_structured_2.msh");
  Mesh mesh = meshgenerator::StructuredMeshGenerator().generate(grid);
  gmsh.write(mesh);

  Log::info() << grid.spec() << std::endl;

  Grid newgrid( grid.spec() );
  Log::info() << newgrid.spec() << std::endl;

  Log::info() << "original: " << grid.uid() << std::endl;
  Log::info() << "fromspec: " << newgrid.uid() << std::endl;
  BOOST_CHECK_EQUAL( grid, newgrid );
}

BOOST_AUTO_TEST_CASE( test_structured_3 )
{
  StructuredGrid grid( "O32" );
  BOOST_CHECK( grid );

  Log::info() << grid.spec() << std::endl;

  Grid newgrid( grid.spec() );
  Log::info() << newgrid.spec() << std::endl;

  Log::info() << "original: " << grid.uid() << std::endl;
  Log::info() << "fromspec: " << newgrid.uid() << std::endl;
  BOOST_CHECK_EQUAL( grid, newgrid );
  BOOST_CHECK_EQUAL( grid.name(), "O32" );
  BOOST_CHECK_EQUAL( newgrid.name(), "O32" );
}


BOOST_AUTO_TEST_CASE( test_iterator )
{
  Grid grid("O4");

  for( atlas::PointXY xy : grid.xy() ) {
    Log::debug() << xy << '\n';
  }

  for( atlas::PointLonLat ll : grid.lonlat() ) {
    Log::debug() << ll << '\n';
  }
  Log::debug() << std::flush;
}


} // namespace test
} // namespace atlas
