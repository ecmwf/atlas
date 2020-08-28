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
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace interpolation {
namespace method {

void KNearestNeighboursBase::buildPointSearchTree( Mesh& meshSource ) {
    eckit::TraceTimer<Atlas> tim( "KNearestNeighboursBase::buildPointSearchTree()" );


    array::ArrayView<double, 2> lonlat = array::make_view<double, 2>( meshSource.nodes().lonlat() );

    static bool fastBuildKDTrees = eckit::Resource<bool>( "$ATLAS_FAST_BUILD_KDTREES", true );

    if ( fastBuildKDTrees ) {
        pTree_.reserve( lonlat.shape( 0 ) );
    }
    for ( idx_t ip = 0; ip < lonlat.shape( 0 ); ++ip ) {
        pTree_.insert( PointLonLat( lonlat( ip, LON ), lonlat( ip, LAT ) ), ip );
    }
    pTree_.build();


    // generate 3D point coordinates
    mesh::actions::BuildXYZField( "xyz" )( meshSource );
}

void KNearestNeighboursBase::buildPointSearchTree( const functionspace::Points& points ) {
    eckit::TraceTimer<Atlas> tim( "KNearestNeighboursBase::buildPointSearchTree()" );

    static bool fastBuildKDTrees = eckit::Resource<bool>( "$ATLAS_FAST_BUILD_KDTREES", true );

    if ( fastBuildKDTrees ) {
        pTree_.reserve( points.size() );
    }
    size_t ip{0};
    for ( auto p : points.iterate().lonlat() ) {
        pTree_.insert( p, ip++ );
    }
    pTree_.build();
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
