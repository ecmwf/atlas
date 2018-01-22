/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/grid/Grid.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_ll_meshgen_one_part" )
{
  Grid g( "L5" );
  Mesh m = meshgenerator::StructuredMeshGenerator().generate(g);
  output::Gmsh("L5.msh").write(m);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    return atlas::test::run( argc, argv );
}

