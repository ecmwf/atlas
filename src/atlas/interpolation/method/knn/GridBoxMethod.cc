/*
 * (C) Copyright 2013 ECMWF.
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
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"


namespace atlas {
namespace interpolation {
namespace method {


MethodBuilder<GridBoxMethod> __builder( "grid-box-average" );


void give_up( const std::forward_list<size_t>& failures ) {
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


GridBoxMethod::GridBoxMethod( const Method::Config& config ) : KNearestNeighboursBase( config ) {
    config.get( "matrix_free", matrixFree_ = false );
    config.get( "fail_early", failEarly_ = true );
}


GridBoxMethod::~GridBoxMethod() = default;


void GridBoxMethod::print( std::ostream& out ) const {
    out << "GridBoxMethod[]";
}


void GridBoxMethod::setup( const FunctionSpace& /*source*/, const FunctionSpace& /*target*/ ) {
    ATLAS_NOTIMPLEMENTED;
}


bool GridBoxMethod::intersect( size_t i, const util::GridBox& box, const PointIndex3::NodeList& closest,
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


void GridBoxMethod::setup( const Grid& source, const Grid& target ) {
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
    ATLAS_ASSERT( pTree_ != nullptr );

    sourceBoxes_ = util::GridBoxes( source );
    targetBoxes_ = util::GridBoxes( target );

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
            if ( ++progress ) {
                Log::info() << "GridBoxMethod: " << *pTree_ << std::endl;
            }

            if ( intersect( i, targetBoxes_.at( i ), pTree_->findInSphere( p, searchRadius_ ), triplets ) ) {
                std::copy( triplets.begin(), triplets.end(), std::back_inserter( allTriplets ) );
            }

            ++i;
        }

        if ( !failures_.empty() ) {
            give_up( failures_ );
        }
    }

    {
        ATLAS_TRACE( "GridBoxMethod::setup: build interpolant matrix" );
        Matrix A( targetBoxes_.size(), sourceBoxes_.size(), allTriplets );
        matrix_.swap( A );
    }
}


void GridBoxMethod::execute( const FieldSet& source, FieldSet& target ) const {
    if ( matrixFree_ ) {
        ATLAS_TRACE( "GridBoxMethod::execute()" );


        // ensure setup()
        functionspace::Points tgt = target_;
        ATLAS_ASSERT( tgt );

        ATLAS_ASSERT( pTree_ != nullptr );
        ATLAS_ASSERT( searchRadius_ > 0. );
        ATLAS_ASSERT( !sourceBoxes_.empty() );
        ATLAS_ASSERT( !targetBoxes_.empty() );


        // set arrays
        ATLAS_ASSERT( source.size() == 1 && source.field( 0 ).rank() == 1 );
        ATLAS_ASSERT( target.size() == 1 && target.field( 0 ).rank() == 1 );

        auto xarray = atlas::array::make_view<double, 1>( source.field( 0 ) );
        auto yarray = atlas::array::make_view<double, 1>( target.field( 0 ) );
        ATLAS_ASSERT( xarray.size() == idx_t( sourceBoxes_.size() ) );
        ATLAS_ASSERT( yarray.size() == idx_t( targetBoxes_.size() ) );

        yarray.assign( 0. );
        failures_.clear();


        // interpolate
        eckit::ProgressTimer progress( "Intersecting", targetBoxes_.size(), "grid box", double( 5. ) );

        std::vector<Triplet> triplets;
        size_t i = 0;
        for ( auto p : tgt.iterate().xyz() ) {
            if ( ++progress ) {
                Log::info() << "GridBoxMethod: " << *pTree_ << std::endl;
            }

            if ( intersect( i, targetBoxes_.at( i ), pTree_->findInSphere( p, searchRadius_ ), triplets ) ) {
                auto& y = yarray[i /*t.col()*/];
                for ( auto& t : triplets ) {
                    y += xarray[t.col()] * t.value();
                }
            }

            ++i;
        }

        if ( !failures_.empty() ) {
            give_up( failures_ );
        }

        return;
    }

    Method::execute( source, target );
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
