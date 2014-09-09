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

#include "tests/TestMeshes.hpp"
#include "atlas/mpl/MPL.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/actions/BuildPeriodicBoundaries.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/actions/BuildEdges.hpp"
#include "atlas/actions/BuildDualMesh.hpp"

BOOST_AUTO_TEST_CASE( test_read_write )
{
  using namespace atlas;

  MPL::init();
  int nlat = 5;
  int lon[5] = {10, 12, 14, 16, 16};
  Mesh::Ptr mesh = test::generate_mesh(nlat, lon);
  Gmsh().write(*mesh,"mesh.msh");

  BOOST_REQUIRE_NO_THROW( mesh = Mesh::Ptr( Gmsh::read( "mesh.msh" ) ) );


  MPL::finalize();
}
