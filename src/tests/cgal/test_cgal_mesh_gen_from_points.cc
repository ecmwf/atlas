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
#include "atlas/Tesselation.h"
#include "atlas/io/Gmsh.h"
#include "atlas/mpi/mpi.h"

//------------------------------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::io;

//------------------------------------------------------------------------------------------------------

#define NLATS 64
#define NLONG 128

//------------------------------------------------------------------------------------------------------

int main()
{
	atlas::mpi::init();

    Mesh::Ptr mesh( new Mesh() );

    Tesselation::generate_lonlat_points( *mesh, NLATS, NLONG );

    Tesselation::tesselate(*mesh);

	Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

	atlas::mpi::finalize();

    return 0;
}
