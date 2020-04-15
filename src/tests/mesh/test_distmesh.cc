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
#include <sstream>

#include "eckit/types/FloatCompare.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/TestMeshes.h"

using namespace atlas;
using namespace atlas::output;
using namespace atlas::meshgenerator;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

double dual_volume( Mesh& mesh ) {
    mesh::Nodes& nodes = mesh.nodes();
    mesh::IsGhostNode is_ghost_node( nodes );
    idx_t nb_nodes                           = nodes.size();
    array::ArrayView<double, 1> dual_volumes = array::make_view<double, 1>( nodes.field( "dual_volumes" ) );
    double area                              = 0;

    for ( idx_t node = 0; node < nb_nodes; ++node ) {
        if ( !is_ghost_node( node ) ) {
            area += dual_volumes( node );
        }
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) { mpi::comm().allReduceInPlace( area, eckit::mpi::sum() ); }

    return area;
}

//-----------------------------------------------------------------------------

CASE( "test_distribute_t63" ) {
    // Every task builds full mesh
    //  meshgenerator::StructuredMeshGenerator generate( util::Config
    //      ("nb_parts",1)
    //      ("part",0) );
    StructuredMeshGenerator generate( util::Config( "partitioner", "equal_regions" ) );

    // long lon[] = {4,6,8,8,8};
    // test::TestGrid grid(5,lon);

    //  GG grid(120,60);
    Grid grid( "N16" );

    Mesh m( generate( grid ) );

    //  actions::distribute_mesh(m);

    mesh::actions::build_parallel_fields( m );
    mesh::actions::build_periodic_boundaries( m );
    mesh::actions::build_halo( m, 1 );
    // mesh::actions::renumber_nodes_glb_idx(m.nodes());
    mesh::actions::build_edges( m );

    Gmsh( "dist.msh", util::Config( "ghost", true ) ).write( m );


    mesh::actions::build_edges_parallel_fields( m );
    mesh::actions::build_median_dual_mesh( m );

    double computed_dual_volume = test::dual_volume( m );
    EXPECT( eckit::types::is_approximately_equal( computed_dual_volume, 360. * 180., 0.0001 ) );
    double difference = 360. * 180. - computed_dual_volume;
    if ( difference > 1e-8 ) {
        std::cout << "difference = " << difference << std::endl;
    }

    Gmsh( "dist.msh" ).write( m );

    mesh::actions::write_load_balance_report( m, "load_balance.dat" );

    Mesh& mesh1 = m;
    EXPECT( mesh1.nodes().size() == m.nodes().size() );

    const array::ArrayView<int, 1> part  = array::make_view<int, 1>( m.nodes().partition() );
    const array::ArrayView<int, 1> ghost = array::make_view<int, 1>( m.nodes().ghost() );
    const array::ArrayView<int, 1> flags = array::make_view<int, 1>( m.nodes().flags() );

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
