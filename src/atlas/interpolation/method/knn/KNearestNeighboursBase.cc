/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "eckit/config/Resource.h"
#include "eckit/log/TraceTimer.h"

#include "atlas/array.h"
#include "atlas/interpolation/method/knn/KNearestNeighboursBase.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace interpolation {
namespace method {

void KNearestNeighboursBase::buildPointSearchTree( Mesh& meshSource ) {
    using namespace atlas;
    eckit::TraceTimer<Atlas> tim( "atlas::interpolation::method::KNearestNeighboursBase::setup()" );

    // generate 3D point coordinates
    mesh::actions::BuildXYZField( "xyz" )( meshSource );
    array::ArrayView<double, 2> coords = array::make_view<double, 2>( meshSource.nodes().field( "xyz" ) );

    // build point-search tree
    pTree_.reset( new PointIndex3 );

    static bool fastBuildKDTrees = eckit::Resource<bool>( "$ATLAS_FAST_BUILD_KDTREES", true );

    if ( fastBuildKDTrees ) {
        std::vector<PointIndex3::Value> pidx;
        pidx.reserve( meshSource.nodes().size() );
        for ( idx_t ip = 0; ip < meshSource.nodes().size(); ++ip ) {
            PointIndex3::Point p{coords( ip, 0 ), coords( ip, 1 ), coords( ip, 2 )};
            pidx.emplace_back( p, ip );
        }
        pTree_->build( pidx.begin(), pidx.end() );
    }
    else {
        for ( idx_t ip = 0; ip < meshSource.nodes().size(); ++ip ) {
            PointIndex3::Point p{coords( ip, 0 ), coords( ip, 1 ), coords( ip, 2 )};
            pTree_->insert( PointIndex3::Value( p, ip ) );
        }
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
