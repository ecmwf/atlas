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

#include "atlas/atlas_config.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/generators/Delaunay.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/util/parallel/mpi/mpi.h"

//------------------------------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::util::io;
using namespace atlas::mesh::generators;

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
    mesh::Mesh::Ptr mesh( generate(*grid) );

    Gmsh gmsh;
    gmsh.options.set<std::string>("nodes","xyz");
    gmsh.write(*mesh, std::string("earth.msh") );

    eckit::mpi::finalize();

    return 0;
}
