/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/interpolation/method/knn/KNearestNeighbours.h"

#include "eckit/log/Plural.h"
#include "eckit/log/Timer.h"

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/interpolation/method/MethodFactory.h"
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

MethodBuilder<KNearestNeighbours> __builder( "k-nearest-neighbours" );

}  // namespace

KNearestNeighbours::KNearestNeighbours( const Method::Config& config ) : KNearestNeighboursBase( config ) {
    k_ = 1;
    config.get( "k-nearest-neighbours", k_ );
    ATLAS_ASSERT( k_ );
}

void KNearestNeighbours::do_setup( const Grid& source, const Grid& target ) {
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

void KNearestNeighbours::do_setup( const FunctionSpace& source, const FunctionSpace& target ) {
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
    weights_triplets.reserve( out_npts * k_ );
    {
        Trace timer( Here(), "atlas::interpolation::method::NearestNeighbour::do_setup()" );

        std::vector<double> weights;

        for ( size_t ip = 0; ip < out_npts; ++ip ) {
            if ( ip && ( ip % 1000 == 0 ) ) {
                double rate = ip / timer.elapsed();
                Log::debug() << eckit::BigNum( ip ) << " (at " << rate << " points/s)..." << std::endl;
            }

            // find the closest input points to the output point
            PointIndex3::Point p{coords( ip, (size_t)0 ), coords( ip, (size_t)1 ), coords( ip, (size_t)2 )};
            PointIndex3::NodeList nn = pTree_->kNearestNeighbours( p, k_ );

            // calculate weights (individual and total, to normalise) using distance
            // squared
            const size_t npts = nn.size();
            ATLAS_ASSERT( npts );
            weights.resize( npts, 0 );

            double sum = 0;
            for ( size_t j = 0; j < npts; ++j ) {
                PointIndex3::Point np = nn[j].point();
                const double d2       = eckit::geometry::Point3::distance2( p, np );

                weights[j] = 1. / ( 1. + d2 );
                sum += weights[j];
            }
            ATLAS_ASSERT( sum > 0 );

            // insert weights into the matrix
            for ( size_t j = 0; j < npts; ++j ) {
                size_t jp = nn[j].payload();
                ATLAS_ASSERT( jp < inp_npts );
                weights_triplets.emplace_back( ip, jp, weights[j] / sum );
            }
        }
    }

    // fill sparse matrix and return
    Matrix A( out_npts, inp_npts, weights_triplets );
    matrix_.swap( A );
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
