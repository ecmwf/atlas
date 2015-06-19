/*
 * (C) Copyright 1996-2014 ECMWF.
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

#include "atlas/atlas_config.h"

#include "atlas/Mesh.h"
#include "atlas/Grid.h"
#include "atlas/meshgen/Delaunay.h"
#include "atlas/io/Gmsh.h"
#include "atlas/mpi/mpi.h"

//------------------------------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::io;
using namespace atlas::meshgen;

//------------------------------------------------------------------------------------------------------

#define NLATS 64
#define NLONG 128

//------------------------------------------------------------------------------------------------------

int main()
{
    eckit::mpi::init();

    Grid::Ptr grid( Grid::create( "ll.32x11") );

    // Build a mesh from grid
    Delaunay generate;
    Mesh::Ptr mesh( generate(*grid) );

    Gmsh gmsh;
    gmsh.options.set<std::string>("nodes","xyz");
    gmsh.write(*mesh, std::string("earth.msh") );

    eckit::mpi::finalize();

    return 0;
}
