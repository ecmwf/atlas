/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestLL
#include "ecbuild/boost_test_framework.h"

#include "atlas/atlas.h"
#include "eckit/mpi/Comm.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/output/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/grid/lonlat/RegularLonLat.h"

#include "tests/AtlasFixture.h"


namespace atlas {
namespace test {

struct AtlasFixture {
    AtlasFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_ll_meshgen_one_part )
{
  grid::lonlat::RegularLonLat g(5);
  mesh::Mesh m;
  mesh::generators::Structured().generate(g,m);
  output::Gmsh("L5.msh").write(m);
}

} // namespace test
} // namespace atlas


