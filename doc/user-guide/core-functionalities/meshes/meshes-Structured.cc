#include "atlas/atlas.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

using atlas::atlas_init;
using atlas::atlas_finalize;
using atlas::grid::Grid;
using atlas::mesh::Mesh;
using atlas::output::Gmsh;
using atlas::util::Config;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    atlas::meshgenerator::StructuredMeshGenerator meshgenerator;

    Grid grid( "O32" );
    Mesh::Ptr mesh( meshgenerator.generate(grid) );

    Gmsh gmsh_2d("mesh2d.msh");
    Gmsh gmsh_3d("mesh3d.msh", Config("coordinates", "xyz") );

    gmsh_2d.write(*mesh);
    gmsh_3d.write(*mesh);

    atlas_finalize();

    return 0;
}
