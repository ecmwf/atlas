/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/Output.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/TestMeshes.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_gmsh_output_1" ) {
    Mesh mesh = test::generate_mesh( Grid( "N32" ) );
    output::Gmsh gmsh( "test_gmsh_output_1.msh" );
    gmsh.write( mesh );
}

CASE( "test_gmsh_output_2" ) {
    Mesh mesh = test::generate_mesh( Grid( "N32" ) );
    atlas::output::GmshFileStream file( "test_gmsh_output_2.msh", "w" );
    Log::warning() << "TODO: Not yet implemented!!! ATLAS-254" << std::endl;
    // output::Gmsh gmsh( file );
    // gmsh.write( mesh );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
