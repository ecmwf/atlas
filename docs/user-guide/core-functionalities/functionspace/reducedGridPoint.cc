#include "atlas/atlas.h"
#include "atlas/grids/grids.h"
#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/io/Gmsh.h"
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
    //double  zdist;// zlon, zlat;
    int jnode = 0;

    // Calculate scalar function
    Field::Ptr scalarField1(fs_rgp->createGlobalField<double>("scalar1"));

    Field::Ptr lonField(fs_rgp->createGlobalField<double>("lon"));
    Field::Ptr latField(fs_rgp->createGlobalField<double>("lat"));
    Field::Ptr distField(fs_rgp->createGlobalField<double>("dist"));

    ArrayView <double,1> scalar1(*scalarField1);
    ArrayView <double,1> zlon(*lonField);
    ArrayView <double,1> zlat(*latField);
    ArrayView <double,1> zdist(*distField);

    // Temporary
    if (eckit::mpi::rank() == 0)
    {
        for (int jlat = 0; jlat < reducedGrid->nlat(); ++jlat)
        {
            for (int jlon =0; jlon < reducedGrid->nlon(jlat); ++jlon)
            {
                zlat(jnode) = reducedGrid->lat(jlat);
                zlat(jnode) = zlat(jnode) * deg2rad;

                zlon(jnode) = reducedGrid->lon(jlat, jlon);
                zlon(jnode) = zlon(jnode) * deg2rad;

                zdist(jnode) = 2.0 * sqrt((cos(zlat(jnode)) * sin((zlon(jnode)-zlonc)/2)) *
                              (cos(zlat(jnode)) * sin((zlon(jnode)-zlonc)/2)) +
                               sin((zlat(jnode)-zlatc)/2) * sin((zlat(jnode)-zlatc)/2));

                scalar1(jnode) = 0.0;
                if (zdist(jnode) < zrad)
                {
                    scalar1(jnode) = 0.5 * (1. + cos(rpi*zdist(jnode)/zrad));
                }
                jnode = jnode+1;
            }
        }
    }

    Field::Ptr localField(fs_rgp->createField<double>("local"));
    fs_rgp->scatter(*scalarField1, *localField);

    // Write mesh and field in gmsh format for visualization
    io::Gmsh gmsh;
    gmsh.options.set("info", true);
    gmsh.write(*mesh, "mesh.msh");
    gmsh.write(*localField, *fs_rgp, "scalar1.msh");
    /* .... */

    if (eckit::mpi::rank() == 0)
    {
        cout << "grid points             = "
             << reducedGrid->npts()   << endl;
        cout << "globalField grid points = "
             << lonField->size()      << endl;
    }

    atlas_finalize();
    return 0;
}
