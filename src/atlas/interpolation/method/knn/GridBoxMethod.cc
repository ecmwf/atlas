/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/knn/GridBoxMethod.h"

#include <algorithm>
#include <vector>

#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"


namespace atlas {
namespace interpolation {
namespace method {


GridBoxMethod::GridBoxMethod( const Method::Config& config ) : KNearestNeighboursBase( config ) {
    config.get( "matrix_free", matrixFree_ = false );
    config.get( "fail_early", failEarly_ = true );
}


GridBoxMethod::~GridBoxMethod() = default;


void GridBoxMethod::print( std::ostream& out ) const {
    out << "GridBoxMethod[]";
}


void GridBoxMethod::do_setup( const FunctionSpace& /*source*/, const FunctionSpace& /*target*/ ) {
    ATLAS_NOTIMPLEMENTED;
}

bool GridBoxMethod::intersect( size_t i, const GridBox& box, const util::IndexKDTree::ValueList& closest,
                               std::vector<eckit::linalg::Triplet>& triplets ) const {
    ASSERT( !closest.empty() );

    triplets.clear();
    triplets.reserve( closest.size() );

    double area = box.area();
    ASSERT( area > 0. );

    double sumSmallAreas = 0.;
    for ( auto& c : closest ) {
        auto j        = c.payload();
        auto smallBox = sourceBoxes_.at( j );

        if ( box.intersects( smallBox ) ) {
            double smallArea = smallBox.area();
            ASSERT( smallArea > 0. );

            triplets.emplace_back( i, j, smallArea / area );
            sumSmallAreas += smallArea;

            if ( eckit::types::is_approximately_equal( area, sumSmallAreas, 1. /*m^2*/ ) ) {
                return true;
            }
        }
    }

    if ( failEarly_ ) {
        Log::error() << "Failed to intersect grid box " << i << ", " << box << std::endl;
        throw_Exception( "Failed to intersect grid box" );
    }

    failures_.push_front( i );
    triplets.clear();
    return false;
}


void GridBoxMethod::do_setup( const Grid& source, const Grid& target ) {
    ATLAS_TRACE( "GridBoxMethod::setup()" );

    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }

    ATLAS_ASSERT( source );
    ATLAS_ASSERT( target );

    functionspace::Points src( source );
    functionspace::Points tgt( target );
    ATLAS_ASSERT( src );
    ATLAS_ASSERT( tgt );
    source_ = src;
    target_ = tgt;

    buildPointSearchTree( src );

    sourceBoxes_ = GridBoxes( source );
    targetBoxes_ = GridBoxes( target );

    searchRadius_ = sourceBoxes_.getLongestGridBoxDiagonal() + targetBoxes_.getLongestGridBoxDiagonal();
    failures_.clear();

    if ( matrixFree_ ) {
        Matrix A;
        matrix_.swap( A );
        return;
    }

    std::vector<Triplet> allTriplets;

    {
        ATLAS_TRACE( "GridBoxMethod::setup: intersecting grid boxes" );

        eckit::ProgressTimer progress( "Intersecting", targetBoxes_.size(), "grid box", double( 5. ) );

        std::vector<Triplet> triplets;
        size_t i = 0;
        for ( auto p : tgt.iterate().xyz() ) {
            ++progress;

            if ( intersect( i, targetBoxes_.at( i ), pTree_.closestPointsWithinRadius( p, searchRadius_ ),
                            triplets ) ) {
                std::copy( triplets.begin(), triplets.end(), std::back_inserter( allTriplets ) );
            }

            ++i;
        }

        if ( !failures_.empty() ) {
            giveUp( failures_ );
        }
    }

    {
        ATLAS_TRACE( "GridBoxMethod::setup: build interpolant matrix" );
        Matrix A( targetBoxes_.size(), sourceBoxes_.size(), allTriplets );
        matrix_.swap( A );
    }
}


void GridBoxMethod::giveUp( const std::forward_list<size_t>& failures ) {
    Log::warning() << "Failed to intersect grid boxes: ";

    size_t count = 0;
    auto sep     = "";
    for ( const auto& f : failures ) {
        if ( count++ < 10 ) {
            Log::warning() << sep << f;
            sep = ", ";
        }
    }
    Log::warning() << "... (" << eckit::Plural( count, "total failure" ) << std::endl;

    throw_Exception( "Failed to intersect grid boxes" );
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
