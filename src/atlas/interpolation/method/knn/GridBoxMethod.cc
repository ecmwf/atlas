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

#include <forward_list>
#include <vector>

#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/functionspace/Points.h"
#include "atlas/grid.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/GridBox.h"


namespace atlas {
namespace interpolation {
namespace method {


MethodBuilder<GridBoxMethod> __builder( "grid-box-average" );


void GridBoxMethod::setup( const Grid& source, const Grid& target ) {
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }

    ATLAS_ASSERT_MSG( !source.projection() && !target.projection(),
                      "GridBoxMethod: rotated/projected grids not supported" );

    sourceGrid_ = source;
    targetGrid_ = target;

    setup( functionspace::Points( source ), functionspace::Points( target ) );
}


GridBoxMethod::~GridBoxMethod() = default;


void GridBoxMethod::print( std::ostream& out ) const {
    out << "GridBoxMethod[]";
}


void GridBoxMethod::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }

    ATLAS_ASSERT( sourceGrid_ );
    ATLAS_ASSERT( targetGrid_ );
    source_ = source;
    target_ = target;

    functionspace::Points src = source;
    functionspace::Points tgt = target;
    ATLAS_ASSERT( src );
    ATLAS_ASSERT( tgt );

    auto src_npts = size_t( src.size() );
    auto tgt_npts = size_t( tgt.size() );
    Log::debug() << "GridBoxMethod: intersect " << eckit::BigNum( tgt_npts ) << " from "
                 << eckit::Plural( tgt_npts, "grid box" ) << std::endl;


    // build point-search tree
    {
        Trace timer( Here(), "GridBoxMethod: build point search tree" );
        buildPointSearchTree( src );
        ATLAS_ASSERT( pTree_ != nullptr );
    }


    // helper structures
    size_t nbFailures = 0;
    std::forward_list<size_t> failures;

    std::vector<Triplet> weights_triplets;
    std::vector<Triplet> triplets;


    // set grid boxes
    const util::GridBoxes sourceBoxes( sourceGrid_ );
    const util::GridBoxes targetBoxes( targetGrid_ );
    const auto R = sourceBoxes.getLongestGridBoxDiagonal() + targetBoxes.getLongestGridBoxDiagonal();


    // intersect grid boxes
    {
        eckit::ProgressTimer progress( "Projecting", tgt_npts, "point", double( 5. ), Log::info() );
        Trace timer( Here(), "GridBoxMethod: intersect grid boxes" );

        size_t i = 0;
        for ( auto p : tgt.iterate().xyz() ) {
            if ( ++progress ) {
                Log::info() << "GridBoxMethod: " << *pTree_ << std::endl;
            }

            // lookup
            auto closest = pTree_->findInSphere( p, R );
            ASSERT( !closest.empty() );


            // calculate grid box intersections
            triplets.clear();
            triplets.reserve( closest.size() );

            auto& box   = targetBoxes.at( i );
            double area = box.area();
            ASSERT( area > 0. );

            double sumSmallAreas = 0.;
            bool areaMatch       = false;
            for ( auto& c : closest ) {
                auto j        = c.payload();
                auto smallBox = sourceBoxes.at( j );

                if ( box.intersects( smallBox ) ) {
                    double smallArea = smallBox.area();
                    ASSERT( smallArea > 0. );

                    triplets.emplace_back( i, j, smallArea / area );
                    sumSmallAreas += smallArea;

                    if ( ( areaMatch = eckit::types::is_approximately_equal( area, sumSmallAreas, 1. /*m^2*/ ) ) ) {
                        break;
                    }
                }
            }


            // insert the interpolant weights into the global (sparse) interpolant matrix
            if ( areaMatch ) {
                std::copy( triplets.begin(), triplets.end(), std::back_inserter( weights_triplets ) );
            }
            else {
                ++nbFailures;
                failures.push_front( i );
            }


            ++i;
        }
    }


    // statistics
    Log::debug() << "Intersected " << eckit::Plural( weights_triplets.size(), "grid box" ) << std::endl;

    if ( nbFailures > 0 ) {
        Log::warning() << "Failed to intersect points: ";
        size_t count = 0;
        auto sep     = "";
        for ( const auto& f : failures ) {
            Log::warning() << sep << f;
            sep = ", ";
            if ( ++count > 10 ) {
                Log::warning() << "...";
                break;
            }
        }
        Log::warning() << std::endl;
    }


    // fill sparse matrix
    {
        Trace timer( Here(), "GridBoxMethod: build sparse matrix" );
        Matrix A( tgt_npts, src_npts, weights_triplets );
        matrix_.swap( A );
    }
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
