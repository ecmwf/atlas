#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/MeshGen.hpp"

//------------------------------------------------------------------------------------------------------

using namespace atlas;

//------------------------------------------------------------------------------------------------------

#define NLATS 256
#define NLONG 256

//------------------------------------------------------------------------------------------------------

int main()
{
    std::vector< Point3 >* pts = atlas::MeshGen::generate_latlon_points(NLATS, NLONG);

    std::cout << "generated " << pts->size() << " points" << std::endl;

    Mesh* mesh = atlas::MeshGen::generate_from_points( *pts );

    atlas::Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

    delete pts;
    delete mesh;

    return 0;
}
