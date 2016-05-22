#include "eckit/mpi/mpi.h"
#include "eckit/config/Resource.h"
#include "atlas/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid/grids.h"
#include "atlas/field/Field.h"
#include "atlas/array/ArrayView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/NodeColumns.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::field;
using namespace atlas::array;
using namespace atlas::grid::global;
using namespace atlas::mesh;


int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Generate global classic reduced Gaussian grid
    string gridID = Resource<string>("--grid", string("N32"));
    SharedPtr<Structured> grid( Structured::create(gridID) );

    // Generate mesh associated to structured grid
    mesh::generators::Structured meshgenerator;
    SharedPtr<Mesh> mesh ( meshgenerator.generate(*grid) );

    // Number of nodes in the mesh
    // (different from number of points on a grid!)
    int nb_nodes = mesh->nodes().size();

    // Number of vertical levels required
    int nb_levels = 10;

    // Generate functionspace associated to mesh
    SharedPtr<functionspace::NodeColumns> fs_nodes(
      new functionspace::NodeColumns(*mesh, mesh::Halo(1)));

    // Note on field generation
    SharedPtr<Field> field_scalar1(
      fs_nodes->createField<double>("scalar1") );
    SharedPtr<Field> field_scalar2(
      fs_nodes->createField<double>("scalar2", nb_levels) );
    SharedPtr<Field> field_vector1(
      fs_nodes->createField<double>("vector1", make_shape(2)) );
    SharedPtr<Field> field_vector2(
      fs_nodes->createField<double>("vector2", nb_levels,
                                               make_shape(2)) );
    SharedPtr<Field> field_tensor1(
      fs_nodes->createField<double>("tensor1", make_shape(2,2)) );
    SharedPtr<Field> field_tensor2(
      fs_nodes->createField<double>("tensor2", nb_levels,
                                               make_shape(2,2)) );
    /* .... */
    // Variables for scalar1 field definition
    const double rpi = 2.0 * asin(1.0);
    const double deg2rad = rpi / 180.;
    const double zlatc = 0.0 * rpi;
    const double zlonc = 1.0 * rpi;
    const double zrad  = 2.0 * rpi / 9.0;
    double  zdist, zlon, zlat;

    // Retrieve lonlat field to calculate scalar1 function
    ArrayView <double,1> scalar1(*field_scalar1);
    ArrayView <double,2> lonlat ( mesh->nodes().lonlat() );
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
    output::Gmsh output_mesh("mesh.msh");
    output_mesh.write(*mesh);
    output::Gmsh output_fields("scalar1.msh");
    output_fields.write(*field_scalar1);

    /* .... */
    // Halo exchange
    fs_nodes->haloExchange(*field_scalar1);
    std::string checksum = fs_nodes->checksum(*field_scalar1);
    Log::info() << checksum << endl;

    // Create a global field
    SharedPtr<Field> field_global(
       fs_nodes->createField("global", *field_scalar1, field::global() ) );
    // Gather operation
    fs_nodes->gather(*field_scalar1, *field_global);

    Log::info() << "local nodes         = "
                << fs_nodes->nb_nodes()  << endl;
    Log::info() << "grid points          = "
                << grid->npts()   << endl;
    Log::info() << "field_global.shape(0) = "
                << field_global->shape(0) << endl;

    // Scatter operation
    fs_nodes->scatter(*field_global, *field_scalar1);

    // Halo exchange and checksum
    fs_nodes->haloExchange(*field_scalar1);
    checksum = fs_nodes->checksum(*field_scalar1);
    Log::info() << checksum << endl;

    // FieldSet checksum
    FieldSet fields;
    fields.add(*field_scalar1);
    fields.add(*field_global);
    checksum = fs_nodes->checksum(fields);
    Log::info() << checksum << endl;
    /* .... */
    // Operations
    size_t N;
    gidx_t gidx_min, gidx_max;
    double min, max, sum, mean, stddev;

    // Minimum and maximum
    fs_nodes->minimum(*field_scalar1, min);
    fs_nodes->maximum(*field_scalar1, max);
    Log::info() << "min: " << min << endl;
    Log::info() << "max: " << max << endl;

    // Minimum and maximum + location
    fs_nodes->minimumAndLocation(*field_scalar1, min, gidx_min);
    fs_nodes->maximumAndLocation(*field_scalar1, max, gidx_max);
    Log::info() << "min: "        << min << ",  "
                << "global_id = " << gidx_min << endl;
    Log::info() << "max: "        << max << ",  "
                << "global_id = " << gidx_max << endl;

    // Summation
    fs_nodes->sum(*field_scalar1, sum, N);
    Log::info() << "sum: "         << sum
                << ", nb_nodes = " << N << endl;

    // Order independent (from partitioning) summation
    fs_nodes->orderIndependentSum(*field_scalar1, sum, N);
    Log::info() << "oi_sum: "      << sum
                << ", nb_nodes = " << N << endl;

    // Average over number of nodes
    fs_nodes->mean(*field_scalar1, mean, N);
    Log::info() << "mean: " << mean << ", nb_nodes = " << N << endl;

    // Average and standard deviation over number of nodes
    fs_nodes->meanAndStandardDeviation(
                *field_scalar1, mean, stddev, N);
    Log::info() << "mean = "         << mean   << ",  "
                << "std_deviation: " << stddev << ",  "
                << "nb_nodes: "      << N      << endl;

    atlas_finalize();

    return 0;
}
