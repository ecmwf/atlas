#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/mesh/actions/GenerateMesh.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/util/io/Gmsh.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace atlas;
using namespace atlas::grids;
using namespace atlas::meshgen;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    string gridID, visualize;
    gridID    = eckit::Resource<string>("--grid"     ,
                                        string("N32"));
    visualize = eckit::Resource<string>("--visualize",
                                        string("2D"));

    ReducedGrid::Ptr reducedGrid(ReducedGrid::create(gridID));

    Mesh::Ptr meshPtr;
    ReducedGridMeshGenerator generate_mesh;
    meshPtr = Mesh::Ptr(generate_mesh(*reducedGrid));
 
    io::Gmsh gmsh;
    gmsh.options.set("info", true);
    if (visualize == "3D")
    {
        actions::BuildXYZField("xyz")(*meshPtr);
        gmsh.options.set("nodes", std::string("xyz"));
    }
    gmsh.write(*meshPtr, "mesh.msh");
    
    atlas_finalize();

    return 0;
}
