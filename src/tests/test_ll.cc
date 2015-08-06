/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestLL
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpi/mpi.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/io/Gmsh.h"
#include "atlas/Mesh.h"
#include "atlas/grids/LonLatGrid.h"


using namespace atlas::io;
using namespace atlas::meshgen;
using namespace atlas::grids;

namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( init ) { eckit::mpi::init(); }

BOOST_AUTO_TEST_CASE( test_ll_meshgen_one_part )
{
  LonLatGrid g(11,LonLatGrid::INCLUDES_POLES);
  Mesh m;
  ReducedGridMeshGenerator().generate(g,m);
  Gmsh().write(m,"lonlat11.msh");
}

BOOST_AUTO_TEST_CASE( finalize ) { eckit::mpi::finalize(); }

} // namespace test
} // namespace atlas


