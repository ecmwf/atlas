/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>

#include "eckit/filesystem/PathName.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Unique.h"

using Topology = atlas::mesh::Nodes::Topology;
using atlas::util::UniqueLonLat;

//------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;
using atlas::util::Config;
using eckit::PathName;

//------------------------------------------------------------------------------
namespace atlas {

struct Node {
    Node( gidx_t gid, int idx ) {
        g = gid;
        i = idx;
    }
    gidx_t g;
    gidx_t i;
    bool operator<( const Node& other ) const { return ( g < other.g ); }
};

//------------------------------------------------------------------------------

void make_nodes_global_index_human_readable( const mesh::actions::BuildHalo& build_halo, mesh::Nodes& nodes,
                                             bool do_all ) {
    ATLAS_TRACE();
    // TODO: ATLAS-14: fix renumbering of EAST periodic boundary points
    // --> Those specific periodic points at the EAST boundary are not checked for
    // uid,
    //     and could receive different gidx for different tasks

    UniqueLonLat compute_uid( nodes );

    // unused // int mypart = mpi::comm().rank();
    int nparts  = mpi::comm().size();
    size_t root = 0;

    array::ArrayView<gidx_t, 1> nodes_glb_idx = array::make_view<gidx_t, 1>( nodes.global_index() );
    // nodes_glb_idx.dump( Log::info() );
    Log::info() << std::endl;
    Log::info() << "min = " << nodes.global_index().metadata().getLong( "min" ) << std::endl;
    Log::info() << "max = " << nodes.global_index().metadata().getLong( "max" ) << std::endl;
    Log::info() << "human_readable = " << nodes.global_index().metadata().getBool( "human_readable" ) << std::endl;
    gidx_t glb_idx_max = 0;

    std::vector<int> points_to_edit;

    if ( do_all ) {
        points_to_edit.resize( nodes_glb_idx.size() );
        for ( idx_t i = 0; i < nodes_glb_idx.size(); ++i ) {
            points_to_edit[i] = i;
        }
    }
    else {
        glb_idx_max = nodes.global_index().metadata().getLong( "max", 0 );
        points_to_edit.resize( build_halo.periodic_points_local_index_.size() );
        for ( size_t i = 0; i < points_to_edit.size(); ++i ) {
            points_to_edit[i] = build_halo.periodic_points_local_index_[i];
        }
    }

    std::vector<gidx_t> glb_idx( points_to_edit.size() );
    for ( size_t i = 0; i < points_to_edit.size(); ++i ) {
        glb_idx[i] = nodes_glb_idx( i );
    }

    ATLAS_DEBUG_VAR( points_to_edit );
    ATLAS_DEBUG_VAR( points_to_edit.size() );

    ATLAS_TRACE( "distributed_sort" );

    /*
 * Sorting following gidx will define global order of
 * gathered fields. Special care needs to be taken for
 * pole edges, as their centroid might coincide with
 * other edges
 */
    int nb_nodes = glb_idx.size();

    // 1) Gather all global indices, together with location

    std::vector<int> recvcounts( mpi::comm().size() );
    std::vector<int> recvdispls( mpi::comm().size() );

    ATLAS_TRACE_MPI( GATHER ) { mpi::comm().gather( nb_nodes, recvcounts, root ); }

    recvdispls[0] = 0;
    for ( int jpart = 1; jpart < nparts; ++jpart )  // start at 1
    {
        recvdispls[jpart] = recvcounts[jpart - 1] + recvdispls[jpart - 1];
    }
    int glb_nb_nodes = std::accumulate( recvcounts.begin(), recvcounts.end(), 0 );

    std::vector<gidx_t> glb_idx_gathered( glb_nb_nodes );
    ATLAS_TRACE_MPI( GATHER ) {
        mpi::comm().gatherv( glb_idx.data(), glb_idx.size(), glb_idx_gathered.data(), recvcounts.data(),
                             recvdispls.data(), root );
    }

    // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
    std::vector<Node> node_sort;
    node_sort.reserve( glb_nb_nodes );
    for ( size_t jnode = 0; jnode < glb_idx_gathered.size(); ++jnode ) {
        node_sort.emplace_back( glb_idx_gathered[jnode], jnode );
    }

    ATLAS_TRACE_SCOPE( "local_sort" ) { std::sort( node_sort.begin(), node_sort.end() ); }

    gidx_t gid = glb_idx_max;
    for ( size_t jnode = 0; jnode < node_sort.size(); ++jnode ) {
        if ( jnode == 0 ) {
            ++gid;
        }
        else if ( node_sort[jnode].g != node_sort[jnode - 1].g ) {
            ++gid;
        }
        int inode               = node_sort[jnode].i;
        glb_idx_gathered[inode] = gid;
    }

    // 3) Scatter renumbered back
    ATLAS_TRACE_MPI( SCATTER ) {
        mpi::comm().scatterv( glb_idx_gathered.data(), recvcounts.data(), recvdispls.data(), glb_idx.data(),
                              glb_idx.size(), root );
    }

    for ( int jnode = 0; jnode < nb_nodes; ++jnode ) {
        nodes_glb_idx( points_to_edit[jnode] ) = glb_idx[jnode];
    }

    // nodes_glb_idx.dump( Log::info() );
    // Log::info() << std::endl;
}

}  // namespace atlas
//-----------------------------------------------------------------------------

class Tool : public AtlasTool {
    int execute( const Args& args ) override;
    std::string briefDescription() override {
        return "Tool to generate a python script that plots the grid-distribution "
               "of a given grid";
    }
    std::string usage() override { return name() + " (--grid=name) [--help]"; }

public:
    Tool( int argc, char** argv );

private:
    std::string key;
    PathName path_in;
    PathName path_out;
};

//-----------------------------------------------------------------------------

Tool::Tool( int argc, char** argv ) : AtlasTool( argc, argv ) {
    add_option( new SimpleOption<std::string>(
        "grid", "Grid unique identifier\n" + indent() + "     Example values: N80, F40, O24, L32" ) );
    add_option( new SimpleOption<long>( "halo", "size of halo" ) );
    add_option( new SimpleOption<bool>( "do-all", "Renumber all points" ) );
}

//-----------------------------------------------------------------------------

int Tool::execute( const Args& args ) {
    Trace t( Here(), "main" );

    key = "";
    args.get( "grid", key );

    std::string path_in_str = "";
    if ( args.get( "grid", path_in_str ) ) {
        path_in = path_in_str;
    }

    StructuredGrid grid;
    if ( key.size() ) {
        try {
            grid = Grid( key );
        }
        catch ( eckit::Exception& ) {
        }
    }
    else {
        Log::error() << "No grid specified." << std::endl;
    }

    if ( !grid ) {
        return failed();
    }

    Log::debug() << "Domain: " << grid.domain() << std::endl;
    Log::debug() << "Periodic: " << grid.periodic() << std::endl;

    MeshGenerator meshgenerator( "structured", Config( "partitioner", "equal_regions" ) );

    size_t halo = args.getLong( "halo", 0 );
    bool do_all = args.getBool( "do-all", false );

    ATLAS_DEBUG_VAR( do_all );

    size_t iterations = 1;
    mpi::comm().barrier();

    for ( size_t j = 0; j < 1; ++j ) {
        ATLAS_TRACE( "outer_iteration" );
        Mesh mesh = meshgenerator.generate( grid );

        atlas::mesh::actions::build_periodic_boundaries( mesh );

        Log::info() << "building halo" << std::endl;
        atlas::mesh::actions::BuildHalo build_halo( mesh );
        build_halo( halo );

        Trace::Barriers set_barrier( true );
        Trace::Tracing set_channel( Log::info() );
        for ( size_t i = 0; i < iterations; ++i ) {
            make_nodes_global_index_human_readable( build_halo, mesh.nodes(), do_all );
        }
    }
    t.stop();
    Log::info() << Trace::report( Config( "indent", 2 )( "decimals", 2 )
                                  // ("depth",5)
    );
    return success();
}

//------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    Tool tool( argc, argv );
    return tool.start();
}
