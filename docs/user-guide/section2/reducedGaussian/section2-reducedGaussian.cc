#include "atlas/atlas.h"
#include "atlas/grids/grids.h"
#include "atlas/io/Gmsh.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/actions/GenerateMesh.h"
#include "atlas/actions/BuildXYZField.h"

using namespace atlas;
using namespace atlas::grids;
using namespace atlas::meshgen;

int main(int argc, char *argv[])
{
    try
    {
        if (argc < 3)
		{
            throw "Mesh not generated! "
                  "The executable needs 2 command line arguments: "
                  "<filename> <dimensions>";
		}
    }
	catch (const char *msg)
	{
        std::cerr << msg << std::endl;
		return EXIT_FAILURE;
	}
    std::string surfdimStr = argv[2];
    int surfdim = atoi(surfdimStr.c_str());
    std::string fileName = argv[1];
    
    atlas_init(argc, argv);

    Grid::Ptr RGgridPtr32(Grid::create("N32"));

    Mesh::Ptr meshPtr;
    ReducedGridMeshGenerator generate_mesh;
    meshPtr = Mesh::Ptr(generate_mesh(grids::rgg::N32()));
 
    io::Gmsh gmsh; 
    gmsh.options.set("info", true);
    if (surfdim == 3)
    {
        actions::BuildXYZField("xyz")(*meshPtr);
        gmsh.options.set("nodes", std::string("xyz"));
    }
    gmsh.write(*meshPtr, fileName);
    
    atlas_finalize();

    return EXIT_SUCCESS;
}



