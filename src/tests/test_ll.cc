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
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/grid/LonLatGrid.h"


using namespace atlas::util::io;
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

BOOST_AUTO_TEST_CASE( test_ll_meshgen_one_part )
{
  LonLatGrid g(11,LonLatGrid::INCLUDES_POLES);
  mesh::Mesh m;
  ReducedGridMeshGenerator().generate(g,m);
  Gmsh().write(m,"lonlat11.msh");
}

} // namespace test
} // namespace atlas


