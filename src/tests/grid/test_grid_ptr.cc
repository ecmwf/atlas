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


#include "atlas/library/atlas.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"


#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"

using Grid       = atlas::grid::Grid;
using Structured = atlas::grid::StructuredGrid;
using Regular    = atlas::grid::RegularGrid;
using Config     = atlas::util::Config;

namespace atlas {
namespace test {

struct GlobalFixture {
    GlobalFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv); }
    ~GlobalFixture() { atlas_finalize(); }
};
BOOST_GLOBAL_FIXTURE( GlobalFixture );

BOOST_AUTO_TEST_CASE( test_from_string_L32 )
{
  Grid grid;
  BOOST_CHECK( not grid );

  grid = Grid("L32");
  BOOST_CHECK( grid );
  BOOST_CHECK( Structured(grid) == true  );
  BOOST_CHECK( Regular(grid)    == true  );

  auto structured = Structured(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 65 );
  BOOST_CHECK_EQUAL( structured.periodic(), true );
  BOOST_CHECK_EQUAL( structured.nx(0), 128 );
  BOOST_CHECK_EQUAL( structured.y().front(), 90. );
  BOOST_CHECK_EQUAL( structured.y().back(), -90. );

  auto regular = Regular(grid);
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

  BOOST_CHECK_EQUAL( Structured(grid), true  );
  BOOST_CHECK_EQUAL( Regular(grid)   , false );

  auto structured = Structured(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 64 );
  BOOST_CHECK_EQUAL( structured.periodic(), true );
  BOOST_CHECK_EQUAL( structured.nx().front(), 20 );
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

  grid = Structured( json_config );
  BOOST_CHECK( grid );
  BOOST_CHECK( Structured(grid) == true  );
  BOOST_CHECK( Regular(grid)    == true  );

  auto structured = Structured(grid);
  BOOST_CHECK_EQUAL( structured.ny(), 9 );
  BOOST_CHECK_EQUAL( structured.periodic(), true );
  BOOST_CHECK_EQUAL( structured.nx(0), 16 );
  BOOST_CHECK_EQUAL( structured.y().front(), 90. );
  BOOST_CHECK_EQUAL( structured.y().back(), -90. );

  auto regular = Regular(grid);
  BOOST_CHECK_EQUAL( regular.ny(), 9 );
  BOOST_CHECK_EQUAL( regular.periodic(), true );
  BOOST_CHECK_EQUAL( regular.nx(), 16 );
  BOOST_CHECK_EQUAL( regular.y().front(), 90. );
  BOOST_CHECK_EQUAL( regular.y().back(), -90. );

  output::Gmsh gmsh("test_grid_ptr.msh");
  eckit::SharedPtr<mesh::Mesh> mesh ( meshgenerator::StructuredMeshGenerator().generate(grid) );
  gmsh.write(*mesh);
}

} // namespace test
} // namespace atlas
