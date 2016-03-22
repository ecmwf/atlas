#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "eckit/config/Resource.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::array;
using namespace atlas::grid::global;
using namespace atlas::field;
using namespace atlas::mesh;


int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Generate global reduced grid
    string gridID = Resource<string>("--grid", string("N32"));
    SharedPtr<Structured> grid (Structured::create(gridID));

    // Number of points in the grid
    int const nb_nodes = grid->npts();

    // Generate functionspace associated to grid
    SharedPtr<functionspace::StructuredColumns>
        fs_rgp(new functionspace::StructuredColumns(*grid));

    /* .... */
    // Variables for scalar1 field definition
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    int jnode = 0;



    // Calculate scalar function
    SharedPtr<Field> field_scalar1(fs_rgp->createField<double>("scalar1"));
    ArrayView <double,1> scalar1(*field_scalar1);

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
            ++jnode;
        }
    }

    // Generate mesh associated to reduced grid
    mesh::generators::Structured meshgenerator;
    SharedPtr<Mesh> mesh (meshgenerator.generate(*grid));

    // Write mesh and field in gmsh format for visualization
    util::io::Gmsh gmsh;
    gmsh.options.set("info", true);
    gmsh.write(*mesh, "mesh.msh");
    gmsh.write(*field_scalar1, "scalar1.msh");
    /* .... */

    atlas_finalize();
    return 0;
}
