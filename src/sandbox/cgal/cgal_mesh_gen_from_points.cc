#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
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
    std::unique_ptr<Mesh> mesh( new Mesh() );

    Tesselation::generate_latlon_points( *mesh, NLATS, NLONG );

    Tesselation::tesselate(*mesh);

    Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

    return 0;
}
