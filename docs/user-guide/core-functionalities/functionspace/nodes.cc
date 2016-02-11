#include "atlas/atlas.h"
#include "atlas/grids/grids.h"
#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/io/Gmsh.h"
#include "atlas/functionspace/Nodes.h"
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
    int nb_nodes = mesh->nodes().size();

    // Number of vertical levels required
    int nb_levels = 10;

    // Generate functionspace associated to mesh
    SharedPtr<functionspace::Nodes> fs_nodes(new functionspace::
                                             Nodes(*mesh, Halo(1)));

    // Note on field generation
    Field::Ptr scalarField1(fs_nodes->createField<double>
                      ("scalar1"));
    Field::Ptr scalarField2(fs_nodes->createField<double>
                      ("scalar2", nb_levels));
    Field::Ptr vectorField1(fs_nodes->createField<double>
                      ("vector1", make_shape(2)));
    Field::Ptr vectorField2(fs_nodes->createField<double>
                      ("vector2", nb_levels, make_shape(2)));
    Field::Ptr tensorField1(fs_nodes->createField<double>
                      ("tensor1", make_shape(2,2)));
    Field::Ptr tensorField2(fs_nodes->createField<double>
                      ("tensor2", nb_levels, make_shape(2,2)));
    /* .... */
    // Variables for scalar1 field definition
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    double  zdist, zlon, zlat;

    // Retrieve lonlat field to calculate scalar1 function
    ArrayView <double,1> scalar1(*scalarField1);
    ArrayView <double,2> lonlat = ArrayView<double,2>(
                                    mesh->nodes().lonlat());
    for (int jnode = 0; jnode < nb_nodes; ++jnode)
    {
        zlon = lonlat(jnode,0) * deg2rad;
        zlat = lonlat(jnode,1) * deg2rad;

        zdist = 2.0 * sqrt((cos(zlat) * sin((zlon-zlonc)/2)) *
                           (cos(zlat) * sin((zlon-zlonc)/2)) +
                        sin((zlat-zlatc)/2) * sin((zlat-zlatc)/2));

        scalar1(jnode) = 0.0;
        if (zdist < zrad)
        {
            scalar1(jnode) = 0.5 * (1. + cos(rpi*zdist/zrad));
        }
    }

    // Write mesh and field in gmsh format for visualization
    io::Gmsh gmsh;
    gmsh.options.set("info", true);
    gmsh.write(*mesh, "mesh.msh");
    gmsh.write(*scalarField1, *fs_nodes, "scalar1.msh");
    /* .... */
    // Halo exchange
    fs_nodes->haloExchange(*scalarField1);
    std::string checksum = fs_nodes->checksum(*scalarField1);
    if (eckit::mpi::rank() == 0) cout << checksum << endl;

    // Create a global field
    Field::Ptr globalField(fs_nodes->createGlobalField(
                               "global", *scalarField1));
    // Gather operation
    fs_nodes->gather(*scalarField1, *globalField);

    if (eckit::mpi::rank() == 0)
    {
        cout << "local nodes         = "
             << fs_nodes->nb_nodes()  << endl;
        cout << "grid points          = "
             << reducedGrid->npts()   << endl;
        cout << "globalField.shape(0) = "
             << globalField->shape(0) << endl;
    }

    // Scatter operation
    fs_nodes->scatter(*globalField, *scalarField1);

    // Halo exchange and checksum
    fs_nodes->haloExchange(*scalarField1);
    checksum = fs_nodes->checksum(*scalarField1);
    if (eckit::mpi::rank() == 0) cout << checksum << endl;

    // FieldSet checksum
    FieldSet fields;
    fields.add(*scalarField1);
    fields.add(*globalField);
    checksum = fs_nodes->checksum(fields);
    if (eckit::mpi::rank() == 0) cout << checksum << endl;
    /* .... */
    // Operations
    size_t N;
    gidx_t gidx_min, gidx_max;
    double min, max, sum, mean, stddev;

    // Minimum and maximum
    fs_nodes->minimum(*scalarField1, min);
    fs_nodes->maximum(*scalarField1, max);
    Log::info() << "min: " << min << endl;
    Log::info() << "max: " << max << endl;

    // Minimum and maximum + location
    fs_nodes->minimumAndLocation(*scalarField1, min, gidx_min);
    fs_nodes->maximumAndLocation(*scalarField1, max, gidx_max);
    Log::info() << "min: "        << min << ",  "
                << "global_id = " << gidx_min << endl;
    Log::info() << "max: "        << max << ",  "
                << "global_id = " << gidx_max << endl;

    // Summation
    fs_nodes->sum(*scalarField1, sum, N);
    Log::info() << "sum: "         << sum
                << ", nb_nodes = " << N << endl;

    // Order independent (from partitioning) summation
    fs_nodes->orderIndependentSum(*scalarField1, sum, N);
    Log::info() << "oi_sum: "      << sum
                << ", nb_nodes = " << N << endl;

    // Average over number of nodes
    fs_nodes->mean(*scalarField1, mean, N);
    Log::info() << "mean: " << mean << ", nb_nodes = " << N << endl;

    // Average and standard deviation over number of nodes
    fs_nodes->meanAndStandardDeviation(
                *scalarField1, mean, stddev, N);
    Log::info() << "mean = "         << mean   << ",  "
                << "std_deviation: " << stddev << ",  "
                << "nb_nodes: "      << N      << endl;

    atlas_finalize();

    return 0;
}
