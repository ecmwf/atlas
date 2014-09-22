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

#include "atlas/io/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/grid/Tesselation.h"

//------------------------------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;

//------------------------------------------------------------------------------------------------------

#define NLATS 256
#define NLONG 256

//------------------------------------------------------------------------------------------------------

int main()
{
    Mesh::Ptr mesh( new Mesh() );

    Tesselation::generate_latlon_points( *mesh, NLATS, NLONG );

    Tesselation::tesselate(*mesh);

    Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

    return 0;
}
