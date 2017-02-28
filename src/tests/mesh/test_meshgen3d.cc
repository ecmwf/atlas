/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestMeshGen3D
#include "ecbuild/boost_test_framework.h"

#include "atlas/internals/atlas_config.h"

#include "atlas/atlas.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/parallel/mpi/mpi.h"


using namespace atlas::output;
using namespace atlas::mesh::generators;
using namespace atlas::grid;

namespace atlas {
namespace test {

struct GlobalFixture {
    GlobalFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv); }
    ~GlobalFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture );

BOOST_AUTO_TEST_CASE( test_create_mesh )
{
	mesh::Mesh::Ptr m ( mesh::Mesh::create() );


  util::Config opts;
  opts.set("3d", true); ///< creates links along date-line
  opts.set("include_pole", true);      ///< triangulate the pole point
  mesh::generators::Structured generate(opts);

  // opts.set("nb_parts",1); // default = 1
  // opts.set("part",    0); // default = 0

    m = generate( Grid("N24") );

	Gmsh("out.msh", util::Config("coordinates","xyz") ).write(*m);
}

} // namespace test
} // namespace atlas
