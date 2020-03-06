#include "atlas/grid/Grid.h"
#include "atlas/library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

using atlas::Grid;
using atlas::Mesh;
using atlas::StructuredMeshGenerator;
using atlas::output::Gmsh;
using atlas::util::Config;

int main( int argc, char* argv[] ) {
    atlas::initialize( argc, argv );

    StructuredMeshGenerator meshgenerator;

    Grid grid( "O32" );
    Mesh mesh = meshgenerator.generate( grid );

    Gmsh gmsh_2d( "mesh2d.msh" );
    Gmsh gmsh_3d( "mesh3d.msh", Config( "coordinates", "xyz" ) );

    gmsh_2d.write( mesh );
    gmsh_3d.write( mesh );

    atlas::finalize();

    return 0;
}
