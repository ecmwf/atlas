/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestGmsh
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/actions/BuildEdges.hpp"
#include "atlas/actions/BuildDualMesh.hpp"
#include "atlas/actions/BuildPeriodicBoundaries.hpp"

BOOST_AUTO_TEST_CASE( test_read_write )
{
  using namespace atlas;

  MPL::init();
  Mesh* mesh;
  BOOST_REQUIRE_NO_THROW( mesh = Gmsh::read( "T47.msh" ) );

  build_periodic_boundaries(*mesh);
  build_edges(*mesh);
  build_dual_mesh(*mesh);

  Gmsh::write(*mesh,"bla.msh");
  delete mesh;
  MPL::finalize();
}
