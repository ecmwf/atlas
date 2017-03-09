#include "atlas/atlas.h"
#include "atlas/grid/Grid.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/StructuredColumns.h"

using atlas::array::ArrayView;
using atlas::atlas_finalize;
using atlas::atlas_init;
using atlas::field::Field;
using atlas::functionspace::StructuredColumns;
using atlas::grid::Grid;
using atlas::mesh::Mesh;
using atlas::output::Gmsh;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Generate global reduced grid
    Grid grid( "N32" );

    // Number of points in the grid
    int const nb_nodes = grid.npts();

    // Generate functionspace associated to grid
    StructuredColumns::Ptr
        fs_rgp(new StructuredColumns(grid));

    /* .... */
    // Variables for scalar1 field definition
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    int jnode = 0;



    // Calculate scalar function
    Field::Ptr field_scalar1(fs_rgp->createField<double>("scalar1"));
    ArrayView <double,1> scalar1(*field_scalar1);

    for (size_t jlat = 0; jlat < fs_rgp->nlat(); ++jlat)
    {
        for (size_t jlon = 0; jlon < fs_rgp->nlon(jlat); ++jlon)
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

    // Write field
    {
      // Generate visualisation mesh associated to grid
      atlas::meshgenerator::StructuredMeshGenerator meshgenerator;
      Mesh::Ptr mesh (meshgenerator.generate(grid));

      Gmsh gmsh("scalar1.msh");
      gmsh.write(*mesh);
      gmsh.write(*field_scalar1);
    }

    atlas_finalize();
    return 0;
}
