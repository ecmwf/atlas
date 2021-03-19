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
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>

#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Unique.h"

using atlas::functionspace::NodeColumns;
using atlas::mesh::Halo;

namespace atlas {
namespace mesh {
namespace actions {

namespace {

void global_bounding_box( const mesh::Nodes& nodes, double min[2], double max[2] ) {
    ATLAS_TRACE();

    auto xy            = array::make_view<double, 2>( nodes.xy() );
    const int nb_nodes = nodes.size();
    min[XX]            = std::numeric_limits<double>::max();
    min[YY]            = std::numeric_limits<double>::max();
    max[XX]            = -std::numeric_limits<double>::max();
    max[YY]            = -std::numeric_limits<double>::max();

    for ( int node = 0; node < nb_nodes; ++node ) {
        min[XX] = std::min( min[XX], xy( node, XX ) );
        min[YY] = std::min( min[YY], xy( node, YY ) );
        max[XX] = std::max( max[XX], xy( node, XX ) );
        max[YY] = std::max( max[YY], xy( node, YY ) );
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) {
        mpi::comm().allReduceInPlace( min, 2, eckit::mpi::min() );
        mpi::comm().allReduceInPlace( max, 2, eckit::mpi::max() );
    }
}

struct Node {
    Node() = default;
    Node( gidx_t gid, idx_t idx ) {
        g = gid;
        i = idx;
    }
    gidx_t g;
    idx_t i;

    bool operator<( const Node& other ) const { return ( g < other.g ); }
};

}  // namespace

array::Array* build_centroids_xy( const mesh::HybridElements&, const Field& xy );

void add_centroid_dual_volume_contribution( Mesh& mesh, array::ArrayView<double, 1>& dual_volumes );
void add_median_dual_volume_contribution_cells( const mesh::HybridElements& cells, const mesh::HybridElements& edges,
                                                const mesh::Nodes& nodes, array::Array& array_dual_volumes );
void add_median_dual_volume_contribution_poles( const mesh::HybridElements& edges, const mesh::Nodes& nodes,
                                                array::Array& array_dual_volumes );
void build_dual_normals( Mesh& mesh );
void make_dual_normals_outward( Mesh& mesh );

void build_median_dual_mesh( Mesh& mesh ) {
    ATLAS_TRACE();
    bool median_dual_mesh = false;
    mesh.metadata().get( "median_dual_mesh", median_dual_mesh );
    if ( median_dual_mesh ) {
        return;
    }
    mesh.metadata().set( "median_dual_mesh", true );

    mesh::Nodes& nodes          = mesh.nodes();
    mesh::HybridElements& edges = mesh.edges();
    Field dual_volumes =
        nodes.add( Field( "dual_volumes", array::make_datatype<double>(), array::make_shape( nodes.size() ) ) );

    if ( !mesh.cells().has_field( "centroids_xy" ) ) {
        mesh.cells().add( Field( "centroids_xy", build_centroids_xy( mesh.cells(), mesh.nodes().xy() ) ) );
    }

    if ( !mesh.edges().has_field( "centroids_xy" ) ) {
        mesh.edges().add( Field( "centroids_xy", build_centroids_xy( mesh.edges(), mesh.nodes().xy() ) ) );
    }

    array::ArrayView<double, 1> v = array::make_view<double, 1>( dual_volumes );
    v.assign( 0. );
    add_median_dual_volume_contribution_cells( mesh.cells(), mesh.edges(), mesh.nodes(), dual_volumes );
    add_median_dual_volume_contribution_poles( mesh.edges(), mesh.nodes(), dual_volumes );

    build_dual_normals( mesh );

    Field skewness = mesh.edges().add(
        Field( "skewness", array::make_datatype<double>(), array::make_shape( mesh.edges().size() ) ) );
    Field alpha =
        mesh.edges().add( Field( "alpha", array::make_datatype<double>(), array::make_shape( mesh.edges().size() ) ) );
    array::make_view<double, 1>( skewness ).assign( 0. );
    array::make_view<double, 1>( alpha ).assign( 0.5 );

    functionspace::NodeColumns nodes_fs( mesh );
    {
        ATLAS_TRACE( "halo-exchange dual_volumes" );
        nodes_fs.haloExchange( nodes.field( "dual_volumes" ) );
    }

    functionspace::EdgeColumns edges_fs( mesh );
    {
        ATLAS_TRACE( "halo-exchange dual_normals" );
        edges_fs.haloExchange( edges.field( "dual_normals" ) );
    }
    make_dual_normals_outward( mesh );
}

array::Array* build_centroids_xy( const mesh::HybridElements& elements, const Field& field_xy ) {
    auto xy                               = array::make_view<double, 2>( field_xy );
    array::Array* array_centroids         = array::Array::create<double>( array::make_shape( elements.size(), 2 ) );
    array::ArrayView<double, 2> centroids = array::make_view<double, 2>( *array_centroids );
    idx_t nb_elems                        = elements.size();
    const mesh::HybridElements::Connectivity& elem_nodes = elements.node_connectivity();
    for ( idx_t e = 0; e < nb_elems; ++e ) {
        centroids( e, XX )               = 0.;
        centroids( e, YY )               = 0.;
        const idx_t nb_nodes_per_elem    = elem_nodes.cols( e );
        const double average_coefficient = 1. / static_cast<double>( nb_nodes_per_elem );
        for ( idx_t n = 0; n < nb_nodes_per_elem; ++n ) {
            centroids( e, XX ) += xy( elem_nodes( e, n ), XX );
            centroids( e, YY ) += xy( elem_nodes( e, n ), YY );
        }
        centroids( e, XX ) *= average_coefficient;
        centroids( e, YY ) *= average_coefficient;
    }
    return array_centroids;
}

void add_median_dual_volume_contribution_cells( const mesh::HybridElements& cells, const mesh::HybridElements& edges,
                                                const mesh::Nodes& nodes, array::Array& array_dual_volumes ) {
    ATLAS_TRACE();

    auto dual_volumes = array::make_view<double, 1>( array_dual_volumes );

    auto xy             = array::make_view<double, 2>( nodes.xy() );
    auto cell_centroids = array::make_view<double, 2>( cells.field( "centroids_xy" ) );
    auto edge_centroids = array::make_view<double, 2>( edges.field( "centroids_xy" ) );
    const mesh::HybridElements::Connectivity& cell_edge_connectivity = cells.edge_connectivity();
    const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
    auto field_flags                                                 = array::make_view<int, 1>( cells.flags() );

    auto patch = [&field_flags]( idx_t e ) {
        using Topology = atlas::mesh::Nodes::Topology;
        return Topology::check( field_flags( e ), Topology::PATCH );
    };

    // special ordering for bit-identical results
    idx_t nb_cells = cells.size();
    std::vector<Node> ordering( nb_cells );
    for ( idx_t jcell = 0; jcell < nb_cells; ++jcell ) {
        ordering[jcell] =
            Node( util::unique_lonlat( cell_centroids( jcell, XX ), cell_centroids( jcell, YY ) ), jcell );
    }
    std::sort( ordering.data(), ordering.data() + nb_cells );

    for ( idx_t jcell = 0; jcell < nb_cells; ++jcell ) {
        idx_t icell = ordering[jcell].i;
        if ( patch( icell ) ) {
            continue;
        }
        double x0 = cell_centroids( icell, XX );
        double y0 = cell_centroids( icell, YY );

        for ( idx_t jedge = 0; jedge < cell_edge_connectivity.cols( icell ); ++jedge ) {
            idx_t iedge = cell_edge_connectivity( icell, jedge );
            double x1   = edge_centroids( iedge, XX );
            double y1   = edge_centroids( iedge, YY );
            for ( idx_t jnode = 0; jnode < 2; ++jnode ) {
                idx_t inode       = edge_node_connectivity( iedge, jnode );
                double x2         = xy( inode, XX );
                double y2         = xy( inode, YY );
                double triag_area = std::abs( x0 * ( y1 - y2 ) + x1 * ( y2 - y0 ) + x2 * ( y0 - y1 ) ) * 0.5;
                dual_volumes( inode ) += triag_area;
            }
        }
    }
}

void add_median_dual_volume_contribution_poles( const mesh::HybridElements& edges, const mesh::Nodes& nodes,
                                                array::Array& array_dual_volumes ) {
    ATLAS_TRACE();

    array::ArrayView<double, 1> dual_volumes = array::make_view<double, 1>( array_dual_volumes );
    auto xy                                  = array::make_view<double, 2>( nodes.xy() );
    auto edge_centroids                      = array::make_view<double, 2>( edges.field( "centroids_xy" ) );
    const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
    const mesh::HybridElements::Connectivity& edge_cell_connectivity = edges.cell_connectivity();

    const idx_t nb_edges = edges.size();
    std::map<idx_t, std::vector<idx_t>> node_to_bdry_edge;
    for ( idx_t jedge = 0; jedge < nb_edges; ++jedge ) {
        if ( edge_cell_connectivity( jedge, 0 ) != edge_cell_connectivity.missing_value() &&
             edge_cell_connectivity( jedge, 1 ) == edge_cell_connectivity.missing_value() ) {
            node_to_bdry_edge[edge_node_connectivity( jedge, 0 )].push_back( jedge );
            node_to_bdry_edge[edge_node_connectivity( jedge, 1 )].push_back( jedge );
        }
    }

    const double tol = 1.e-6;
    double min[2], max[2];
    global_bounding_box( nodes, min, max );

    std::map<idx_t, std::vector<idx_t>>::iterator it;
    for ( it = node_to_bdry_edge.begin(); it != node_to_bdry_edge.end(); ++it ) {
        const idx_t jnode              = ( *it ).first;
        std::vector<idx_t>& bdry_edges = ( *it ).second;
        const double x0                = xy( jnode, XX );
        const double y0                = xy( jnode, YY );
        double x1, y1, y2;
        for ( idx_t jedge = 0; jedge < static_cast<idx_t>( bdry_edges.size() ); ++jedge ) {
            const idx_t iedge = bdry_edges[jedge];
            x1                = edge_centroids( iedge, XX );
            y1                = edge_centroids( iedge, YY );

            y2 = 0.;
            if ( std::abs( y1 - max[YY] ) < tol ) {
                y2 = 90.;
            }
            else if ( std::abs( y1 - min[YY] ) < tol ) {
                y2 = -90.;
            }

            if ( y2 != 0 ) {
                const double quad_area = std::abs( ( x1 - x0 ) * ( y2 - y0 ) );
                dual_volumes( jnode ) += quad_area;
            }
        }
    }
}

void build_dual_normals( Mesh& mesh ) {
    ATLAS_TRACE();

    array::ArrayView<double, 2> elem_centroids = array::make_view<double, 2>( mesh.cells().field( "centroids_xy" ) );

    mesh::Nodes& nodes          = mesh.nodes();
    mesh::HybridElements& edges = mesh.edges();
    const idx_t nb_edges        = edges.size();

    array::ArrayView<double, 2> node_xy = array::make_view<double, 2>( nodes.xy() );
    double min[2], max[2];
    global_bounding_box( nodes, min, max );
    double tol = 1.e-6;

    double xl, yl, xr, yr;
    array::ArrayView<double, 2> edge_centroids = array::make_view<double, 2>( edges.field( "centroids_xy" ) );
    array::ArrayView<double, 2> dual_normals   = array::make_view<double, 2>(
        edges.add( Field( "dual_normals", array::make_datatype<double>(), array::make_shape( nb_edges, 2 ) ) ) );

    const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
    const mesh::HybridElements::Connectivity& edge_cell_connectivity = edges.cell_connectivity();

    std::map<idx_t, std::vector<idx_t>> node_to_bdry_edge;
    for ( idx_t jedge = 0; jedge < nb_edges; ++jedge ) {
        if ( edge_cell_connectivity( jedge, 0 ) != edge_cell_connectivity.missing_value() &&
             edge_cell_connectivity( jedge, 1 ) == edge_cell_connectivity.missing_value() ) {
            node_to_bdry_edge[edge_node_connectivity( jedge, 0 )].push_back( jedge );
            node_to_bdry_edge[edge_node_connectivity( jedge, 1 )].push_back( jedge );
        }
    }

    for ( idx_t edge = 0; edge < nb_edges; ++edge ) {
        if ( edge_cell_connectivity( edge, 0 ) == edge_cell_connectivity.missing_value() ) {
            // this is a pole edge
            // only compute for one node
            for ( idx_t n = 0; n < 2; ++n ) {
                idx_t node                     = edge_node_connectivity( edge, n );
                std::vector<idx_t>& bdry_edges = node_to_bdry_edge[node];
                double x[2];
                idx_t cnt                 = 0;
                const idx_t nb_bdry_edges = static_cast<idx_t>( bdry_edges.size() );
                for ( idx_t jedge = 0; jedge < nb_bdry_edges; ++jedge ) {
                    idx_t bdry_edge = bdry_edges[jedge];
                    if ( std::abs( edge_centroids( bdry_edge, YY ) - max[YY] ) < tol ) {
                        edge_centroids( edge, YY ) = 90.;
                        x[cnt]                     = edge_centroids( bdry_edge, XX );
                        ++cnt;
                    }
                    else if ( std::abs( edge_centroids( bdry_edge, YY ) - min[YY] ) < tol ) {
                        edge_centroids( edge, YY ) = -90.;
                        x[cnt]                     = edge_centroids( bdry_edge, XX );
                        ++cnt;
                    }
                }
                if ( cnt == 2 ) {
                    dual_normals( edge, XX ) = 0;
                    if ( node_xy( node, YY ) < 0. ) {
                        dual_normals( edge, YY ) = -std::abs( x[1] - x[0] );
                    }
                    else if ( node_xy( node, YY ) > 0. ) {
                        dual_normals( edge, YY ) = std::abs( x[1] - x[0] );
                    }

                    // std::cout << "pole dual_normal = " << dual_normals(YY,edge) <<
                    // std::endl;
                    break;
                }
            }
        }
        else {
            idx_t left_elem  = edge_cell_connectivity( edge, 0 );
            idx_t right_elem = edge_cell_connectivity( edge, 1 );
            xl               = elem_centroids( left_elem, XX );
            yl               = elem_centroids( left_elem, YY );
            if ( right_elem == edge_cell_connectivity.missing_value() ) {
                xr = edge_centroids( edge, XX );
                yr = edge_centroids( edge, YY );
                ;
                if ( std::abs( yr - max[YY] ) < tol ) {
                    yr = 90.;
                }
                else if ( std::abs( yr - min[YY] ) < tol ) {
                    yr = -90.;
                }
            }
            else {
                xr = elem_centroids( right_elem, XX );
                yr = elem_centroids( right_elem, YY );
            }

            dual_normals( edge, XX ) = yl - yr;
            dual_normals( edge, YY ) = -xl + xr;
        }
    }
}

void make_dual_normals_outward( Mesh& mesh ) {
    ATLAS_TRACE();
    mesh::Nodes& nodes                  = mesh.nodes();
    array::ArrayView<double, 2> node_xy = array::make_view<double, 2>( nodes.xy() );

    mesh::HybridElements& edges                                      = mesh.edges();
    const mesh::HybridElements::Connectivity& edge_cell_connectivity = edges.cell_connectivity();
    const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
    array::ArrayView<double, 2> dual_normals = array::make_view<double, 2>( edges.field( "dual_normals" ) );
    const idx_t nb_edges                     = edges.size();

    for ( idx_t edge = 0; edge < nb_edges; ++edge ) {
        if ( edge_cell_connectivity( edge, 0 ) != edge_cell_connectivity.missing_value() ) {
            // Make normal point from node 1 to node 2
            const idx_t ip1 = edge_node_connectivity( edge, 0 );
            const idx_t ip2 = edge_node_connectivity( edge, 1 );
            double dx       = node_xy( ip2, XX ) - node_xy( ip1, XX );
            double dy       = node_xy( ip2, YY ) - node_xy( ip1, YY );
            if ( dx * dual_normals( edge, XX ) + dy * dual_normals( edge, YY ) < 0 ) {
                dual_normals( edge, XX ) = -dual_normals( edge, XX );
                dual_normals( edge, YY ) = -dual_normals( edge, YY );
            }
        }
    }
}

void build_brick_dual_mesh( const Grid& grid, Mesh& mesh ) {
    auto g = StructuredGrid( grid );
    if ( g ) {
        if ( mpi::size() != 1 ) {
            throw_Exception( "Cannot build_brick_dual_mesh with more than 1 task", Here() );
        }

        mesh::Nodes& nodes = mesh.nodes();
        auto xy            = array::make_view<double, 2>( nodes.xy() );
        auto dual_volumes  = array::make_view<double, 1>(
            nodes.add( Field( "dual_volumes", array::make_datatype<double>(), array::make_shape( nodes.size() ) ) ) );
        auto gidx = array::make_view<gidx_t, 1>( nodes.global_index() );

        int c = 0;
        int n = 0;
        for ( idx_t jlat = 0; jlat < g.ny(); ++jlat ) {
            double lat  = g.y( jlat );
            double latN = ( jlat == 0 ) ? 90. : 0.5 * ( lat + g.y( jlat - 1 ) );
            double latS = ( jlat == g.ny() - 1 ) ? -90. : 0.5 * ( lat + g.y( jlat + 1 ) );
            double dlat = ( latN - latS );
            double dlon = 360. / static_cast<double>( g.nx( jlat ) );

            for ( idx_t jlon = 0; jlon < g.nx( jlat ); ++jlon ) {
                while ( gidx( c ) != n + 1 ) {
                    c++;
                }
                ATLAS_ASSERT( xy( c, XX ) == g.x( jlon, jlat ) );
                ATLAS_ASSERT( xy( c, YY ) == lat );
                dual_volumes( c ) = dlon * dlat;
                ++n;
            }
        }

        functionspace::NodeColumns nodes_fs( mesh );
        nodes_fs.haloExchange( nodes.field( "dual_volumes" ) );
    }
    else {
        throw_Exception( "Cannot build_brick_dual_mesh with mesh provided grid type", Here() );
    }
}

void build_centroid_dual_mesh( Mesh& ) {
    ATLAS_NOTIMPLEMENTED;
    // This requires code below which has not been ported yet
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_median_dual_mesh( Mesh::Implementation* mesh ) {
    ATLAS_ASSERT( mesh != nullptr, "Cannot access uninitialised atlas_Mesh" );
    Mesh m( mesh );
    build_median_dual_mesh( m );
}

void atlas__build_centroid_dual_mesh( Mesh::Implementation* mesh ) {
    ATLAS_ASSERT( mesh != nullptr, "Cannot access uninitialised atlas_Mesh" );
    Mesh m( mesh );
    build_centroid_dual_mesh( m );
}
// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
