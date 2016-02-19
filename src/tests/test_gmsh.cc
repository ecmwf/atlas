/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/atlas.h"
#include "atlas/util/parallel/mpi/mpi.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/private/Debug.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildDualMesh.h"

namespace atlas {
namespace test {

struct AtlasFixture {
    AtlasFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_read_write )
{
    using namespace atlas;
    using namespace atlas::io;

    // Mesh::Ptr mesh = test::generate_mesh(nlat, lon);
    Mesh::Ptr mesh = test::generate_mesh(grids::rgg::N128());

    Gmsh gmsh;
    gmsh.options.set("ascii",true);
    gmsh.write(*mesh,"mesh.msh");

    BOOST_REQUIRE_NO_THROW( mesh = Mesh::Ptr( Gmsh().read( "mesh.msh" ) ) );
}

} // namespace test
} // namespace atlas
