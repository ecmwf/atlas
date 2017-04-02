#include "atlas/library/Library.h"
#include "atlas/grid/Grid.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/StructuredColumns.h"

using atlas::grid::Grid;
using atlas::mesh::Mesh;
using atlas::meshgenerator::StructuredMeshGenerator;
using atlas::functionspace::StructuredColumns;
using atlas::field::Field;
using atlas::array::ArrayView;
using atlas::array::make_view;
using atlas::output::Gmsh;

int main(int argc, char *argv[])
{
    atlas::Library::instance().initialise(argc, argv);

    // Generate global reduced grid
    Grid grid( "N32" );

    // Generate functionspace associated to grid
    StructuredColumns fs_rgp(grid);

    // Variables for scalar1 field definition
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    int jnode = 0;

    // Calculate scalar function
    Field field_scalar1 = fs_rgp.createField<double>("scalar1");
    auto scalar1 = make_view<double,1>(field_scalar1);

    for (size_t jlat = 0; jlat < fs_rgp.ny(); ++jlat)
    {
        for (size_t jlon = 0; jlon < fs_rgp.nx(jlat); ++jlon)
        {
            double zlat = fs_rgp.y(jlat);
            zlat = zlat * deg2rad;
            double zlon = fs_rgp.x(jlon, jlat);
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
      StructuredMeshGenerator meshgenerator;
      Mesh mesh = meshgenerator.generate(grid);

      Gmsh gmsh("scalar1.msh");
      gmsh.write(mesh);
      gmsh.write(field_scalar1);
    }

    atlas::Library::instance().finalise();
    return 0;
}
