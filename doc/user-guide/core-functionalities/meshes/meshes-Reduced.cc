#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/mesh/actions/GenerateMesh.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/util/io/Gmsh.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::mesh::generators;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    string gridID    = Resource<string>("--grid"     , string("N32"));
    string visualize = Resource<string>("--visualize", string("2D") );

    SharedPtr<ReducedGrid> reducedGrid( ReducedGrid::create(gridID) );

    ReducedGridMeshGenerator meshgenerator;
    SharedPtr<Mesh> mesh( meshgenerator.generate(*reducedGrid) );
 
    util::io::Gmsh gmsh;
    gmsh.options.set("info", true);
    if (visualize == "3D")
    {
        actions::BuildXYZField("xyz")(*mesh);
        gmsh.options.set("nodes", std::string("xyz"));
    }
    gmsh.write(*mesh, "mesh.msh");
    
    atlas_finalize();

    return 0;
}
