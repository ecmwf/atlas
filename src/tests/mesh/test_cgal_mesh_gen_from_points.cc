/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::meshgenerator;
using namespace atlas::output;

//----------------------------------------------------------------------------------------------------------------------

#define NLATS 64
#define NLONG 128

int main( int argc, char** argv ) {
    atlas::Library::instance().initialise( argc, argv );
    Grid grid( "L33x11" );

    // Build a mesh from grid
    MeshGenerator generate( "delaunay" );
    Mesh mesh = generate( grid );

    Gmsh gmsh( "earth.msh", util::Config( "coordinates", "xyz" ) );
    gmsh.write( mesh );

    atlas::Library::instance().finalise();
    return 0;
}
