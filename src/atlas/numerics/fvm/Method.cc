/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include "atlas/numerics/fvm/Method.h"
#include <cmath>
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "eckit/exception/Exceptions.h"

// =======================================================

using namespace atlas::mesh::actions;

namespace atlas {
namespace numerics {
namespace fvm {

namespace {

mesh::Halo get_halo( const eckit::Parametrisation& params ) {
    size_t halo_size( 1 );
    params.get( "halo", halo_size );
    return mesh::Halo( halo_size );
}

size_t get_levels( const eckit::Parametrisation& params ) {
    size_t levels( 0 );
    params.get( "levels", levels );
    return levels;
}

double get_radius( const eckit::Parametrisation& params ) {
    double radius = util::Earth::radiusInMeters();
    params.get( "radius", radius );
    return radius;
}

}  // namespace

Method::Method( Mesh& mesh ) : Method::Method( mesh, util::NoConfig() ) {
    setup();
}

Method::Method( Mesh& mesh, const mesh::Halo& halo ) : Method::Method( mesh, util::Config( "halo", halo.size() ) ) {
    setup();
}

Method::Method( Mesh& mesh, const eckit::Configuration& params ) :
    mesh_( mesh ),
    levels_( get_levels( params ) ),
    halo_( get_halo( params ) ),
    nodes_( mesh.nodes() ),
    edges_( mesh.edges() ),
    radius_( get_radius( params ) ) {
    setup();
}

void Method::setup() {
    util::Config node_columns_config;
    node_columns_config.set( "halo", halo_.size() );
    if ( levels_ ) node_columns_config.set( "levels", levels_ );
    node_columns_ = functionspace::NodeColumns( mesh(), node_columns_config );
    if ( edges_.size() == 0 ) {
        build_edges( mesh() );
        build_pole_edges( mesh() );
        build_edges_parallel_fields( mesh() );
        build_median_dual_mesh( mesh() );
        build_node_to_edge_connectivity( mesh() );

        const size_t nnodes = nodes_.size();

        // Compute sign
        {
            const array::ArrayView<int, 1> is_pole_edge = array::make_view<int, 1>( edges_.field( "is_pole_edge" ) );

            const mesh::Connectivity& node_edge_connectivity           = nodes_.edge_connectivity();
            const mesh::MultiBlockConnectivity& edge_node_connectivity = edges_.node_connectivity();
            if ( !nodes_.has_field( "node2edge_sign" ) ) {
                nodes_.add( Field( "node2edge_sign", array::make_datatype<double>(),
                                   array::make_shape( nnodes, node_edge_connectivity.maxcols() ) ) );
            }
            array::ArrayView<double, 2> node2edge_sign =
                array::make_view<double, 2>( nodes_.field( "node2edge_sign" ) );

            atlas_omp_parallel_for( size_t jnode = 0; jnode < nnodes; ++jnode ) {
                for ( size_t jedge = 0; jedge < node_edge_connectivity.cols( jnode ); ++jedge ) {
                    size_t iedge = node_edge_connectivity( jnode, jedge );
                    size_t ip1   = edge_node_connectivity( iedge, 0 );
                    if ( jnode == ip1 )
                        node2edge_sign( jnode, jedge ) = 1.;
                    else {
                        node2edge_sign( jnode, jedge ) = -1.;
                        if ( is_pole_edge( iedge ) ) node2edge_sign( jnode, jedge ) = 1.;
                    }
                }
            }
        }

        // Metrics
        if ( 0 ) {
            const size_t nedges                          = edges_.size();
            const array::ArrayView<double, 2> lonlat_deg = array::make_view<double, 2>( nodes_.lonlat() );
            array::ArrayView<double, 1> dual_volumes = array::make_view<double, 1>( nodes_.field( "dual_volumes" ) );
            array::ArrayView<double, 2> dual_normals = array::make_view<double, 2>( edges_.field( "dual_normals" ) );

            const double deg2rad = M_PI / 180.;
            atlas_omp_parallel_for( size_t jnode = 0; jnode < nnodes; ++jnode ) {
                double y  = lonlat_deg( jnode, LAT ) * deg2rad;
                double hx = radius_ * std::cos( y );
                double hy = radius_;
                double G  = hx * hy;
                dual_volumes( jnode ) *= std::pow( deg2rad, 2 ) * G;
            }

            atlas_omp_parallel_for( size_t jedge = 0; jedge < nedges; ++jedge ) {
                dual_normals( jedge, LON ) *= deg2rad;
                dual_normals( jedge, LAT ) *= deg2rad;
            }
        }
    }
    edge_columns_ = functionspace::EdgeColumns( mesh() );
}

// ------------------------------------------------------------------------------------------
extern "C" {
Method* atlas__numerics__fvm__Method__new( Mesh::Implementation* mesh, const eckit::Configuration* params ) {
    Method* method( 0 );
    ATLAS_ERROR_HANDLING( ASSERT( mesh ); Mesh m( mesh ); method = new Method( m, *params ); );
    return method;
}

const functionspace::detail::NodeColumns* atlas__numerics__fvm__Method__functionspace_nodes( Method* This ) {
    ATLAS_ERROR_HANDLING(
        ASSERT( This ); return dynamic_cast<const functionspace::detail::NodeColumns*>( This->node_columns().get() ); );
    return 0;
}

const functionspace::detail::EdgeColumns* atlas__numerics__fvm__Method__functionspace_edges( Method* This ) {
    ATLAS_ERROR_HANDLING(
        ASSERT( This ); return dynamic_cast<const functionspace::detail::EdgeColumns*>( This->edge_columns().get() ); );
    return 0;
}
}
// ------------------------------------------------------------------------------------------

}  // namespace fvm
}  // namespace numerics
}  // namespace atlas
