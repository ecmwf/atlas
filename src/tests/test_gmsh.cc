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
#include "ecbuild/boost_test_framework.h"

#include "tests/TestMeshes.h"
#include "atlas/mpi/mpi.h"
#include "atlas/io/Gmsh.h"
#include "atlas/util/Debug.h"
#include "atlas/Mesh.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildDualMesh.h"

BOOST_AUTO_TEST_CASE( test_read_write )
{
	using namespace atlas;
	using namespace atlas::io;

	eckit::mpi::init();
	int nlat = 5;
	int lon[5] = {10, 12, 14, 16, 16};

	//	Mesh::Ptr mesh = test::generate_mesh(nlat, lon);
	Mesh::Ptr mesh = test::generate_mesh(grids::rgg::N128());

	Gmsh gmsh;
	gmsh.options.set("ascii",true);
	gmsh.write(*mesh,"mesh.msh");

	BOOST_REQUIRE_NO_THROW( mesh = Mesh::Ptr( Gmsh::read( "mesh.msh" ) ) );
	eckit::mpi::finalize();
}
