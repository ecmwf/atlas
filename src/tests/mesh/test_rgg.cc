/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "eckit/types/FloatCompare.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/grid.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Metadata.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace grid {
namespace spacing {
namespace gaussian {
void compute_gaussian_quadrature_npole_equator( const size_t N, double lats[], double weights[] );
}
}  // namespace spacing
}  // namespace grid
}  // namespace atlas

#define DISABLE if ( 0 )
#define ENABLE if ( 1 )

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

static ReducedGaussianGrid debug_grid() {
    return {6, 10, 18, 22, 22, 22, 22, 18, 10, 6};
}

static StructuredGrid minimal_grid( int N, long lon[] ) {
    std::vector<long> nx( 2 * N );
    for ( long j = 0; j < N; ++j ) {
        nx[j]                 = lon[j];
        nx[nx.size() - 1 - j] = nx[j];
    }
    return ReducedGaussianGrid( nx );
}

double compute_lonlat_area( Mesh& mesh ) {
    mesh::Nodes& nodes                 = mesh.nodes();
    mesh::Elements& quads              = mesh.cells().elements( 0 );
    mesh::Elements& triags             = mesh.cells().elements( 1 );
    array::ArrayView<double, 2> lonlat = array::make_view<double, 2>( nodes.lonlat() );

    const mesh::BlockConnectivity& quad_nodes  = quads.node_connectivity();
    const mesh::BlockConnectivity& triag_nodes = triags.node_connectivity();

    double area = 0;
    for ( idx_t e = 0; e < quads.size(); ++e ) {
        idx_t n0  = quad_nodes( e, 0 );
        idx_t n1  = quad_nodes( e, 1 );
        idx_t n2  = quad_nodes( e, 2 );
        idx_t n3  = quad_nodes( e, 3 );
        double x0 = lonlat( n0, LON ), x1 = lonlat( n1, LON ), x2 = lonlat( n2, LON ), x3 = lonlat( n3, LON );
        double y0 = lonlat( n0, LAT ), y1 = lonlat( n1, LAT ), y2 = lonlat( n2, LAT ), y3 = lonlat( n3, LAT );
        area += std::abs( x0 * ( y1 - y2 ) + x1 * ( y2 - y0 ) + x2 * ( y0 - y1 ) ) * 0.5;
        area += std::abs( x2 * ( y3 - y0 ) + x3 * ( y0 - y2 ) + x0 * ( y2 - y3 ) ) * 0.5;
    }
    for ( idx_t e = 0; e < triags.size(); ++e ) {
        idx_t n0  = triag_nodes( e, 0 );
        idx_t n1  = triag_nodes( e, 1 );
        idx_t n2  = triag_nodes( e, 2 );
        double x0 = lonlat( n0, LON ), x1 = lonlat( n1, LON ), x2 = lonlat( n2, LON );
        double y0 = lonlat( n0, LAT ), y1 = lonlat( n1, LAT ), y2 = lonlat( n2, LAT );
        area += std::abs( x0 * ( y1 - y2 ) + x1 * ( y2 - y0 ) + x2 * ( y0 - y1 ) ) * 0.5;
    }
    return area;
}

//-----------------------------------------------------------------------------

CASE( "test_eq_caps" ) {
    std::vector<int> n_regions;
    std::vector<double> s_cap;

    grid::detail::partitioner::eq_caps( 6, n_regions, s_cap );
    EXPECT( n_regions.size() == 3 );
    EXPECT( n_regions[0] == 1 );
    EXPECT( n_regions[1] == 4 );
    EXPECT( n_regions[2] == 1 );

    grid::detail::partitioner::eq_caps( 10, n_regions, s_cap );
    EXPECT( n_regions.size() == 4 );
    EXPECT( n_regions[0] == 1 );
    EXPECT( n_regions[1] == 4 );
    EXPECT( n_regions[2] == 4 );
    EXPECT( n_regions[3] == 1 );
}

CASE( "test_partitioner" ) {
    Grid g( "S4x2" );

    // 12 partitions
    {
        grid::detail::partitioner::EqualRegionsPartitioner partitioner( 12 );
        EXPECT( partitioner.nb_bands() == 4 );
        EXPECT( partitioner.nb_regions( 0 ) == 1 );
        EXPECT( partitioner.nb_regions( 1 ) == 5 );
        EXPECT( partitioner.nb_regions( 2 ) == 5 );
        EXPECT( partitioner.nb_regions( 3 ) == 1 );
    }

    // 24 partitions
    {
        grid::detail::partitioner::EqualRegionsPartitioner partitioner( 24 );
        EXPECT( partitioner.nb_bands() == 5 );
        EXPECT( partitioner.nb_regions( 0 ) == 1 );
        EXPECT( partitioner.nb_regions( 1 ) == 6 );
        EXPECT( partitioner.nb_regions( 2 ) == 10 );
        EXPECT( partitioner.nb_regions( 3 ) == 6 );
        EXPECT( partitioner.nb_regions( 4 ) == 1 );
    }

    // 48 partitions
    {
        grid::detail::partitioner::EqualRegionsPartitioner partitioner( 48 );
        EXPECT( partitioner.nb_bands() == 7 );
        EXPECT( partitioner.nb_regions( 0 ) == 1 );
        EXPECT( partitioner.nb_regions( 1 ) == 6 );
        EXPECT( partitioner.nb_regions( 2 ) == 11 );
        EXPECT( partitioner.nb_regions( 3 ) == 12 );
        EXPECT( partitioner.nb_regions( 4 ) == 11 );
        EXPECT( partitioner.nb_regions( 5 ) == 6 );
        EXPECT( partitioner.nb_regions( 6 ) == 1 );
    }

    // 96 partitions
    {
        grid::detail::partitioner::EqualRegionsPartitioner partitioner( 96 );
        EXPECT( partitioner.nb_bands() == 10 );
        EXPECT( partitioner.nb_regions( 0 ) == 1 );
        EXPECT( partitioner.nb_regions( 1 ) == 6 );
        EXPECT( partitioner.nb_regions( 2 ) == 11 );
        EXPECT( partitioner.nb_regions( 3 ) == 14 );
        EXPECT( partitioner.nb_regions( 4 ) == 16 );
        EXPECT( partitioner.nb_regions( 5 ) == 16 );
        EXPECT( partitioner.nb_regions( 6 ) == 14 );
        EXPECT( partitioner.nb_regions( 7 ) == 11 );
        EXPECT( partitioner.nb_regions( 8 ) == 6 );
        EXPECT( partitioner.nb_regions( 9 ) == 1 );
    }
}

CASE( "test_gaussian_latitudes" ) {
    std::vector<double> factory_latitudes;
    std::vector<double> computed_latitudes;
    std::vector<double> computed_weights;

    size_t size_test_N = 19;

    size_t test_N[] = {16,  24,  32,  48,  64,  80,   96,   128,  160,  200,  256, 320,
                       400, 512, 576, 640, 800, 1024, 1280, 1600, 2000, 4000, 8000};

    for ( size_t i = 0; i < size_test_N; ++i ) {
        size_t N = test_N[i];
        Log::info() << "Testing gaussian latitude " << N << std::endl;
        factory_latitudes.resize( N );
        computed_latitudes.resize( N );
        computed_weights.resize( N );
        // grid::gaussian::latitudes::gaussian_latitudes_npole_equator (N,
        // factory_latitudes.data());
        // grid::gaussian::latitudes::compute_gaussian_quadrature_npole_equator(N,
        // computed_latitudes.data(), computed_weights.data());
        grid::spacing::gaussian::gaussian_latitudes_npole_equator( N, factory_latitudes.data() );
        grid::spacing::gaussian::compute_gaussian_quadrature_npole_equator( N, computed_latitudes.data(),
                                                                            computed_weights.data() );
        double wsum = 0;
        for ( size_t i = 0; i < N; ++i ) {
            EXPECT( eckit::types::is_approximately_equal( computed_latitudes[i], factory_latitudes[i], 0.0000001 ) );
            wsum += computed_weights[i];
        }
        EXPECT( eckit::types::is_approximately_equal( wsum * 2., 1., 0.0000001 ) );
    }
}

CASE( "test_rgg_meshgen_one_part" ) {
    Mesh m;
    util::Config default_opts;
    default_opts.set( "nb_parts", 1 );
    default_opts.set( "part", 0 );
    //  generate.options.set("nb_parts",1);
    //  generate.options.set("part",    0);
    DISABLE {  // This is all valid for meshes generated with MINIMAL NB TRIAGS
        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", true )( "include_pole", false ) );
            m = generate( atlas::test::debug_grid() );
            EXPECT( m.nodes().size() == 156 );
            EXPECT( m.cells().elements( 0 ).size() == 134 );
            EXPECT( m.cells().elements( 1 ).size() == 32 );
            //    EXPECT( m.nodes().metadata().get<size_t>("nb_owned") ==    156 );
            //    EXPECT( m.function_space("quads"
            //    ).metadata().get<size_t>("nb_owned") ==    134 );
            //    EXPECT(
            //    m.function_space("triags").metadata().get<size_t>("nb_owned") ==
            //    32 );
        }

        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", false )( "include_pole", false ) );
            m = generate( atlas::test::debug_grid() );
            EXPECT( m.nodes().size() == 166 );
            EXPECT( m.cells().elements( 0 ).size() == 134 );
            EXPECT( m.cells().elements( 1 ).size() == 32 );
            //    EXPECT( m.nodes().metadata().get<size_t>("nb_owned") ==    166 );
            //    EXPECT( m.function_space("quads"
            //    ).metadata().get<size_t>("nb_owned") ==    134 );
            //    EXPECT(
            //    m.function_space("triags").metadata().get<size_t>("nb_owned") ==
            //    32 );
        }

        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", true )( "include_pole", true ) );
            m = generate( atlas::test::debug_grid() );
            EXPECT( m.nodes().size() == 158 );
            EXPECT( m.cells().elements( 0 ).size() == 134 );
            EXPECT( m.cells().elements( 1 ).size() == 44 );
            //    EXPECT( m.nodes().metadata().get<size_t>("nb_owned") ==    158 );
            //    EXPECT( m.function_space("quads"
            //    ).metadata().get<size_t>("nb_owned") ==    134 );
            //    EXPECT(
            //    m.function_space("triags").metadata().get<size_t>("nb_owned") ==
            //    44 );
        }

        Mesh mesh;

        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", false )( "include_pole", false ) );
            int nlat   = 2;
            long lon[] = {4, 6};
            mesh       = generate( test::minimal_grid( nlat, lon ) );
            EXPECT( mesh.nodes().size() == 24 );
            EXPECT( mesh.cells().elements( 0 ).size() == 14 );
            EXPECT( mesh.cells().elements( 1 ).size() == 4 );

            double max_lat = test::minimal_grid( nlat, lon ).y().front();
            EXPECT( eckit::types::is_approximately_equal( test::compute_lonlat_area( mesh ), 2. * M_PI * 2. * max_lat,
                                                          1e-8 ) );
            output::Gmsh( "minimal2.msh" ).write( mesh );
        }
        // 3 latitudes
        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", false )( "include_pole", false ) );
            int nlat   = 3;
            long lon[] = {4, 6, 8};
            mesh       = generate( test::minimal_grid( nlat, lon ) );
            EXPECT( mesh.nodes().size() == 42 );
            EXPECT( mesh.cells().elements( 0 ).size() == 28 );
            EXPECT( mesh.cells().elements( 1 ).size() == 8 );
            output::Gmsh( "minimal3.msh" ).write( mesh );
        }
        // 4 latitudes
        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", false )( "include_pole", false ) );
            int nlat   = 4;
            long lon[] = {4, 6, 8, 10};
            mesh       = generate( test::minimal_grid( nlat, lon ) );
            EXPECT( mesh.nodes().size() == 64 );
            EXPECT( mesh.cells().elements( 0 ).size() == 46 );
            EXPECT( mesh.cells().elements( 1 ).size() == 12 );
            output::Gmsh( "minimal4.msh" ).write( mesh );
        }
        // 5 latitudes WIP
        ENABLE {
            StructuredMeshGenerator generate( default_opts( "3d", false )( "include_pole", false ) );
            int nlat   = 5;
            long lon[] = {6, 10, 18, 22, 22};
            mesh       = generate( test::minimal_grid( nlat, lon ) );
            EXPECT( mesh.nodes().size() == 166 );
            EXPECT( mesh.cells().elements( 0 ).size() == 134 );
            EXPECT( mesh.cells().elements( 1 ).size() == 32 );
            output::Gmsh( "minimal5.msh" ).write( mesh );
        }
    }
}

CASE( "test_rgg_meshgen_many_parts" ) {
    EXPECT( grid::detail::partitioner::PartitionerFactory::has( "equal_regions" ) );
    int nb_parts = 20;
    //  Alternative grid for debugging
    //  int nlat=10;
    //  long lon[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
    //  test::MinimalGrid grid(nlat,lon);
    StructuredGrid grid = Grid( "N32" );
    // RegularGrid grid(128,64);

    /*
std::cout << grid.spec() << std::endl;
for (int jlat=0;jlat<2*nlat; jlat++) {
      std::cout << grid.lon(jlat,0) << ", ... , " << "," << grid.lon(jlat,9) <<
grid.lon(jlat,10) << std::endl;
}
ASSERT(0);
*/

    double max_lat    = grid.y().front();
    double check_area = 360. * 2. * max_lat;
    double area       = 0;
    int nodes[]  = {313, 332, 336, 338, 334, 337, 348, 359, 360, 361, 360, 360, 359, 370, 321, 334, 338, 335, 348, 315};
    int quads[]  = {243, 290, 293, 294, 291, 293, 310, 320, 321, 321, 320, 321, 320, 331, 278, 291, 294, 293, 305, 245};
    int triags[] = {42, 13, 12, 13, 12, 14, 0, 1, 0, 1, 1, 0, 1, 0, 14, 12, 13, 11, 14, 42};
    int nb_owned = 0;

    std::vector<int> all_owned( grid.size(), -1 );

    for ( int p = 0; p < nb_parts; ++p ) {
        ATLAS_DEBUG_VAR( p );

        StructuredMeshGenerator generate( util::Config( "partitioner", "equal_regions" )( "nb_parts", nb_parts )(
            "part", p )( "include_pole", false )( "3d", false ) );
        ATLAS_DEBUG_HERE();

        Mesh m = generate( grid );
        ATLAS_DEBUG_HERE();
        m.metadata().set( "part", p );
        Log::info() << "generated grid " << p << std::endl;
        array::ArrayView<int, 1> part    = array::make_view<int, 1>( m.nodes().partition() );
        array::ArrayView<gidx_t, 1> gidx = array::make_view<gidx_t, 1>( m.nodes().global_index() );

        area += test::compute_lonlat_area( m );
        ATLAS_DEBUG_HERE();

        DISABLE {  // This is all valid for meshes generated with MINIMAL NB TRIAGS
            if ( nb_parts == 20 ) {
                EXPECT( m.nodes().size() == nodes[p] );
                EXPECT( m.cells().elements( 0 ).size() == quads[p] );
                EXPECT( m.cells().elements( 1 ).size() == triags[p] );
            }
        }
        ATLAS_DEBUG_HERE();

        output::Gmsh( "T63.msh" ).write( m );

        mesh::Nodes& nodes = m.nodes();
        idx_t nb_nodes     = nodes.size();

        // Test if all nodes are connected
        {
            std::vector<int> node_elem_connections( nb_nodes, 0 );

            const mesh::HybridElements::Connectivity& cell_node_connectivity = m.cells().node_connectivity();
            for ( idx_t jelem = 0; jelem < static_cast<idx_t>( m.cells().size() ); ++jelem ) {
                for ( idx_t jnode = 0; jnode < cell_node_connectivity.cols( jelem ); ++jnode ) {
                    node_elem_connections[cell_node_connectivity( jelem, jnode )]++;
                }
            }
            for ( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
                if ( node_elem_connections[jnode] == 0 ) {
                    std::stringstream ss;
                    ss << "part " << p << ": node_gid " << gidx( jnode ) << " is not connected to any element.";
                    DISABLE { Log::error() << ss.str() << std::endl; }
                }
            }
        }

        // Test if all nodes are owned
        for ( idx_t n = 0; n < nb_nodes; ++n ) {
            if ( gidx( n ) <= grid.size() ) {
                if ( part( n ) == p ) {
                    ++nb_owned;
                    if ( all_owned[gidx( n ) - 1] != -1 ) {
                        std::cout << "node " << gidx( n ) << " already visited by " << all_owned[gidx( n ) - 1]
                                  << std::endl;
                    }
                    EXPECT( all_owned[gidx( n ) - 1] == -1 );
                    all_owned[gidx( n ) - 1] = part( n );
                }
            }
        }
    }

    for ( size_t gid = 1; gid <= all_owned.size(); ++gid ) {
        if ( all_owned[gid - 1] == -1 ) {
            Log::error() << "node " << gid << " is not owned by anyone" << std::endl;
        }
    }
    EXPECT( nb_owned == grid.size() );

    EXPECT( eckit::types::is_approximately_equal( area, check_area, 1e-10 ) );
}

CASE( "test_meshgen_ghost_at_end" ) {
    ATLAS_DEBUG_HERE();

    Grid grid( "O8" );

    atlas::util::Config cfg;
    cfg.set( "part", 1 );
    cfg.set( "nb_parts", 8 );
    StructuredMeshGenerator meshgenerator( cfg );
    Mesh mesh        = meshgenerator.generate( grid );
    const auto part  = array::make_view<int, 1>( mesh.nodes().partition() );
    const auto ghost = array::make_view<int, 1>( mesh.nodes().ghost() );
    const auto flags = array::make_view<int, 1>( mesh.nodes().flags() );

    Log::info() << "partition = [ ";
    for ( idx_t jnode = 0; jnode < part.size(); ++jnode ) {
        Log::info() << part( jnode ) << " ";
    }
    Log::info() << "]" << std::endl;

    Log::info() << "ghost     = [ ";
    for ( idx_t jnode = 0; jnode < part.size(); ++jnode ) {
        Log::info() << ghost( jnode ) << " ";
    }
    Log::info() << "]" << std::endl;

    Log::info() << "flags     = [ ";
    for ( idx_t jnode = 0; jnode < part.size(); ++jnode ) {
        Log::info() << mesh::Nodes::Topology::check( flags( jnode ), mesh::Nodes::Topology::GHOST ) << " ";
        EXPECT( mesh::Nodes::Topology::check( flags( jnode ), mesh::Nodes::Topology::GHOST ) == ghost( jnode ) );
    }
    Log::info() << "]" << std::endl;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
