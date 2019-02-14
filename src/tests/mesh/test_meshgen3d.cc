/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::output;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_create_mesh" ) {
    Mesh m;

    util::Config opts;
    opts.set( "3d", true );            ///< creates links along date-line
    opts.set( "include_pole", true );  ///< triangulate the pole point
    StructuredMeshGenerator generate( opts );

    // opts.set("nb_parts",1); // default = 1
    // opts.set("part",    0); // default = 0

    m = generate( Grid( "N24" ) );

    Grid grid = m.grid();
    std::cout << grid.spec() << std::endl;

    Gmsh( "out.msh", util::Config( "coordinates", "xyz" ) ).write( m );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
