/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "eckit/log/Plural.h"

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/knn/NearestNeighbour.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<NearestNeighbour> __builder( "nearest-neighbour" );

}  // namespace

void NearestNeighbour::do_setup( const Grid& source, const Grid& target ) {
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }
    auto functionspace = []( const Grid& grid ) -> FunctionSpace {
        Mesh mesh;
        if ( StructuredGrid( grid ) ) {
            mesh = MeshGenerator( "structured", util::Config( "3d", true ) ).generate( grid );
        }
        else {
            mesh = MeshGenerator( "delaunay" ).generate( grid );
        }
        return functionspace::NodeColumns( mesh );
    };

    do_setup( functionspace( source ), functionspace( target ) );
}

void NearestNeighbour::do_setup( const FunctionSpace& source, const FunctionSpace& target ) {
    source_                        = source;
    target_                        = target;
    functionspace::NodeColumns src = source;
    functionspace::NodeColumns tgt = target;
    ATLAS_ASSERT( src );
    ATLAS_ASSERT( tgt );

    Mesh meshSource = src.mesh();
    Mesh meshTarget = tgt.mesh();

    // build point-search tree
    buildPointSearchTree( meshSource );
    ATLAS_ASSERT( pTree_ != nullptr );

    // generate 3D point coordinates
    mesh::actions::BuildXYZField( "xyz" )( meshTarget );
    array::ArrayView<double, 2> coords = array::make_view<double, 2>( meshTarget.nodes().field( "xyz" ) );

    size_t inp_npts = meshSource.nodes().size();
    size_t out_npts = meshTarget.nodes().size();

    // fill the sparse matrix
    std::vector<Triplet> weights_triplets;
    weights_triplets.reserve( out_npts );
    {
        Trace timer( Here(), "atlas::interpolation::method::NearestNeighbour::do_setup()" );
        for ( size_t ip = 0; ip < out_npts; ++ip ) {
            if ( ip && ( ip % 1000 == 0 ) ) {
                double rate = ip / timer.elapsed();
                Log::debug() << eckit::BigNum( ip ) << " (at " << rate << " points/s)..." << std::endl;
            }

            // find the closest input point to the output point
            PointIndex3::Point p{coords( ip, (size_t)0 ), coords( ip, (size_t)1 ), coords( ip, (size_t)2 )};
            PointIndex3::NodeInfo nn = pTree_->nearestNeighbour( p );
            size_t jp                = nn.payload();

            // insert the weights into the interpolant matrix
            ATLAS_ASSERT( jp < inp_npts );
            weights_triplets.emplace_back( ip, jp, 1 );
        }
    }

    // fill sparse matrix and return
    Matrix A( out_npts, inp_npts, weights_triplets );
    matrix_.swap( A );
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
