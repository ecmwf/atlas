#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "atlas/field/Field.h"
#include "atlas/util/array/ArrayView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/functionspace/ReducedGridPoint.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace atlas;
using namespace eckit;
using namespace atlas::grids;
using namespace atlas::meshgen;
using namespace atlas::functionspace;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Generate global reduced grid
    string gridID = Resource<string>("--grid", string("N32"));
    ReducedGrid::Ptr reducedGrid(ReducedGrid::create(gridID));

    // Generate mesh associated to reduced grid
    ReducedGridMeshGenerator generate_mesh;
    Mesh::Ptr mesh = Mesh::Ptr(generate_mesh(*reducedGrid));

    // Number of nodes in the mesh
    // (different from number of points on a grid!)
    int const nb_nodes = reducedGrid->npts();

    // Generate functionspace associated to mesh
    SharedPtr<functionspace::ReducedGridPoint>
        fs_rgp(new functionspace::ReducedGridPoint(*reducedGrid));//, Halo(1)));

    /* .... */
    // Variables for scalar1 field definition
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    int jnode = 0;



    // Calculate scalar function
    Field::Ptr scalarField1(fs_rgp->createField<double>("scalar1"));
    ArrayView <double,1> scalar1(*scalarField1);

    for (int jlat = 0; jlat < fs_rgp->nlat(); ++jlat)
    {
        for (int jlon = 0; jlon < fs_rgp->nlon(jlat); ++jlon)
        {
            double zlat = fs_rgp->lat(jlat);
            zlat = zlat * deg2rad;
            double zlon = fs_rgp->lon(jlat, jlon);
            zlon  = zlon * deg2rad;
            double zdist = 2.0 * sqrt((cos(zlat) * sin((zlon-zlonc)/2.)) *
                          (cos(zlat) * sin((zlon-zlonc)/2.)) +
                           sin((zlat-zlatc)/2.) * sin((zlat-zlatc)/2.));

            scalar1(jnode) = 0.0;
            if (zdist < zrad)
            {
                scalar1(jnode) = 0.5 * (1. + cos(rpi*zdist/zrad));
            }
            jnode = jnode+1;
        }
    }

    // Write mesh and field in gmsh format for visualization
    io::Gmsh gmsh;
    gmsh.options.set("info", true);
    gmsh.write(*mesh, "mesh.msh");
    gmsh.write(*scalarField1, *fs_rgp, "scalar1.msh");
    /* .... */

    atlas_finalize();
    return 0;
}
