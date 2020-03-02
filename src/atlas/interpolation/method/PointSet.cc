/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/config/Resource.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/interpolation/method/PointSet.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"

using namespace eckit;

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

PointSet::PointSet( const std::vector<Point>& ipts ) : npts_( ipts.size() ) {
    build( ipts );
}

PointSet::PointSet( Mesh& mesh ) {
    mesh::Nodes& nodes = mesh.nodes();

    npts_ = nodes.size();

    ATLAS_ASSERT( npts_ > 0 );

    ATLAS_ASSERT( nodes.has_field( "xyz" ) );

    array::ArrayView<double, 2> coords = array::make_view<double, 2>( nodes.field( "xyz" ) );
    static bool fastBuildKDTrees       = eckit::Resource<bool>( "$ATLAS_FAST_BUILD_KDTREES", true );

    tree_ = new PointIndex3();

    if ( fastBuildKDTrees ) {
        std::vector<PointIndex3::Value> pidx;
        pidx.reserve( npts_ );

        for ( size_t ip = 0; ip < npts_; ++ip ) {
            pidx.emplace_back(
                PointIndex3::Point( coords( ip, (size_t)0 ), coords( ip, (size_t)1 ), coords( ip, (size_t)2 ) ), ip );
        }
        tree_->build( pidx.begin(), pidx.end() );
    }
    else {
        for ( size_t ip = 0; ip < npts_; ++ip ) {
            tree_->insert( PointIndex3::Value(
                PointIndex3::Point( coords( ip, (size_t)0 ), coords( ip, (size_t)1 ), coords( ip, (size_t)2 ) ), ip ) );
        }
    }
}

size_t PointSet::search_unique( const Point& p, size_t idx, uint32_t n ) {
    PointIndex3::NodeList nearest = tree_->kNearestNeighbours( p, Kn( n ) );

    std::vector<size_t> equals;
    equals.reserve( nearest.size() );

    for ( size_t i = 0; i < nearest.size(); ++i ) {
        Point np    = nearest[i].value().point();
        size_t nidx = nearest[i].value().payload();

        //            std::cout << "      - " << nidx << " " << np << std::endl;

        if ( eckit::geometry::points_equal( p, np ) ) {
            //                std::cout << "      EQUAL !!" << std::endl;
            equals.push_back( nidx );
        }
        else {
            break;
        }
    }

    if ( equals.size() == nearest.size() ) /* need to increase the search to find
                                          all duplicates of this point */
    {
        return this->search_unique( p, idx, ++n );
    }
    else /* stop recursion */
    {
        size_t ret = idx; /* if nothing found return same idx */

        if ( equals.size() >= 1 ) { /* if an equal point was found return the first one */
            ret = equals[0];
        }

        for ( size_t n = 1; n < equals.size(); ++n ) {
            duplicates_[equals[n]] = ret;
        }

        return ret;
    }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
