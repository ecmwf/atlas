/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "atlas/internals/atlas_config.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/generators/Delaunay.h"
#include "atlas/output/Gmsh.h"

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh::generators;
using namespace atlas::output;

//----------------------------------------------------------------------------------------------------------------------

#define NLATS 64
#define NLONG 128

int main()
{
    Grid::Ptr grid( Grid::create( "L32x11") );

    // Build a mesh from grid
    Delaunay generate;
    mesh::Mesh::Ptr mesh( generate(*grid) );

    Gmsh gmsh("earth.msh", util::Config("coordinates","xyz") );
    gmsh.write(*mesh);

    return 0;
}
