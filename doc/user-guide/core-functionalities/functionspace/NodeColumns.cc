#include "atlas/functionspace/NodeColumns.h"
#include "atlas/array/ArrayView.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"

using namespace atlas;
using atlas::gidx_t;
using atlas::idx_t;
using atlas::StructuredGrid;
using atlas::StructuredMeshGenerator;
using atlas::array::make_shape;
using atlas::array::make_view;
using atlas::functionspace::NodeColumns;
using atlas::output::Gmsh;

int main( int argc, char* argv[] ) {
    atlas::initialise( argc, argv );

    // Generate global classic reduced Gaussian grid
    StructuredGrid grid( "N32" );

    // Generate mesh associated to structured grid
    StructuredMeshGenerator meshgenerator;
    Mesh mesh = meshgenerator.generate( grid );

    // Number of nodes in the mesh
    // (different from number of points on a grid!)
    size_t nb_nodes = mesh.nodes().size();

    // Number of vertical levels required
    size_t nb_levels = 10;

    // Generate functionspace associated to mesh
    NodeColumns fs_nodes( mesh, option::levels( nb_levels ) | option::halo( 1 ) );

    // Note on field generation
    Field field_scalar1 = fs_nodes.createField<double>( option::name( "scalar1" ) | option::levels( false ) );
    Field field_scalar2 = fs_nodes.createField<double>( option::name( "scalar2" ) );
    Field field_vector1 =
        fs_nodes.createField<double>( option::name( "vector1" ) | option::levels( false ) | option::variables( 2 ) );
    Field field_vector2 = fs_nodes.createField<double>( option::name( "vector2" ) | option::variables( 2 ) );
    Field field_tensor1 = fs_nodes.createField<double>( option::name( "tensor1" ) | option::levels( false ) |
                                                        option::variables( 2 * 2 ) );
    Field field_tensor2 = fs_nodes.createField<double>( option::name( "tensor2" ) | option::variables( 2 * 2 ) );
    /* .... */
    // Variables for scalar1 field definition
    const double rpi     = 2.0 * asin( 1.0 );
    const double deg2rad = rpi / 180.;
    const double zlatc   = 0.0 * rpi;
    const double zlonc   = 1.0 * rpi;
    const double zrad    = 2.0 * rpi / 9.0;
    double zdist, zlon, zlat;

    // Retrieve lonlat field to calculate scalar1 function
    auto scalar1 = make_view<double, 1>( field_scalar1 );
    auto lonlat  = make_view<double, 2>( mesh.nodes().lonlat() );
    for ( size_t jnode = 0; jnode < nb_nodes; ++jnode ) {
        zlon = lonlat( jnode, size_t( 0 ) ) * deg2rad;
        zlat = lonlat( jnode, size_t( 1 ) ) * deg2rad;

        zdist =
            2.0 * sqrt( ( cos( zlat ) * sin( ( zlon - zlonc ) / 2 ) ) * ( cos( zlat ) * sin( ( zlon - zlonc ) / 2 ) ) +
                        sin( ( zlat - zlatc ) / 2 ) * sin( ( zlat - zlatc ) / 2 ) );

        scalar1( jnode ) = 0.0;
        if ( zdist < zrad ) {
            scalar1( jnode ) = 0.5 * ( 1. + cos( rpi * zdist / zrad ) );
        }
    }

    // Write mesh and field in gmsh format for visualization
    Gmsh gmsh( "scalar1.msh" );
    gmsh.write( mesh );
    gmsh.write( field_scalar1 );

    /* .... */
    // Halo exchange
    fs_nodes.haloExchange( field_scalar1 );
    std::string checksum = fs_nodes.checksum( field_scalar1 );
    Log::info() << checksum << std::endl;

    // Create a global field
    Field field_global = fs_nodes.createField( field_scalar1, option::name( "global" ) | option::global() );
    // Gather operation
    fs_nodes.gather( field_scalar1, field_global );

    Log::info() << "local nodes         = " << fs_nodes.nb_nodes() << std::endl;
    Log::info() << "grid points          = " << grid.size() << std::endl;
    Log::info() << "field_global.shape(0) = " << field_global.shape( 0 ) << std::endl;

    // Scatter operation
    fs_nodes.scatter( field_global, field_scalar1 );

    // Halo exchange and checksum
    fs_nodes.haloExchange( field_scalar1 );
    checksum = fs_nodes.checksum( field_scalar1 );
    Log::info() << field_scalar1.name() << " checksum : " << checksum << std::endl;

    // FieldSet checksum
    FieldSet fields;
    fields.add( field_scalar1 );
    fields.add( field_vector1 );
    checksum = fs_nodes.checksum( fields );
    Log::info() << "FieldSet checksum : " << checksum << std::endl;

    /* .... */

    // Operations
    idx_t N;
    gidx_t gidx_min, gidx_max;
    double min, max, sum, mean, stddev;

    // Minimum and maximum
    fs_nodes.minimum( field_scalar1, min );
    fs_nodes.maximum( field_scalar1, max );
    Log::info() << "min: " << min << std::endl;
    Log::info() << "max: " << max << std::endl;

    // Minimum and maximum + location
    fs_nodes.minimumAndLocation( field_scalar1, min, gidx_min );
    fs_nodes.maximumAndLocation( field_scalar1, max, gidx_max );
    Log::info() << "min: " << min << ",  "
                << "global_id = " << gidx_min << std::endl;
    Log::info() << "max: " << max << ",  "
                << "global_id = " << gidx_max << std::endl;

    // Summation
    fs_nodes.sum( field_scalar1, sum, N );
    Log::info() << "sum: " << sum << ", nb_nodes = " << N << std::endl;

    // Order independent (from partitioning) summation
    fs_nodes.orderIndependentSum( field_scalar1, sum, N );
    Log::info() << "oi_sum: " << sum << ", nb_nodes = " << N << std::endl;

    // Average over number of nodes
    fs_nodes.mean( field_scalar1, mean, N );
    Log::info() << "mean: " << mean << ", nb_nodes = " << N << std::endl;

    // Average and standard deviation over number of nodes
    fs_nodes.meanAndStandardDeviation( field_scalar1, mean, stddev, N );
    Log::info() << "mean = " << mean << ",  "
                << "std_deviation: " << stddev << ",  "
                << "nb_nodes: " << N << std::endl;

    atlas::finalise();
    atlas::mpi::finalize();

    return 0;
}
