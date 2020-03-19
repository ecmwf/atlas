/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/container/KDTree.h"

#include "atlas/library/config.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

/// @brief k-dimensional tree constructable both with 2D (lon,lat) points as with 3D (x,y,z) points
///
/// The implementation is based on eckit::KDTreeMemory with 3D (x,y,z) points.
/// 2D points (lon,lat) are converted when needed to 3D during insertion, and during search.
template <typename PayloadT>
class KDTree {
private:
    struct eckit_KDTreeTraits {
        using Point   = Point3;
        using Payload = PayloadT;
    };
    using eckit_KDTree = eckit::KDTreeMemory<eckit_KDTreeTraits>;

public:
    using Payload  = PayloadT;
    using Node     = typename eckit_KDTree::NodeInfo;
    using NodeList = typename eckit_KDTree::NodeList;

    /// @brief Reserve memory for building the kdtree in one shot (optional, at cost of extra memory)
    void reserve( idx_t size ) { tmp_.reserve( size ); }

    /// @brief Insert spherical point (lon,lat)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void insert( const Point2& p2, const Payload& payload ) {
        Point3 p3;
        util::Earth::convertSphericalToCartesian( p2, p3 );
        insert( p3, payload );
    }

    /// @brief Insert 3D cartesian point (x,y,z)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void insert( const Point3& p3, const Payload& payload ) {
        if ( tmp_.capacity() ) {
            tmp_.emplace_back( p3, payload );
        }
        else {
            tree_.insert( {p3, payload} );
        }
    }

    /// @brief Build the kd-tree in one shot, if memory has been reserved.
    /// This will need to be called before all search functions like kNearestNeighbours().
    void build() {
        if ( tmp_.size() ) {
            tree_.build( tmp_.begin(), tmp_.end() );
            tmp_.clear();
            tmp_.shrink_to_fit();
        }
    }

    /// @brief Find k nearest neighbour given a 2D lonlat point (lon,lat)
    NodeList kNearestNeighbours( const Point2& p2, size_t k ) const {
        Point3 p3;
        util::Earth::convertSphericalToCartesian( p2, p3 );
        return tree_.kNearestNeighbours( p3, k );
    }

    /// @brief Find k nearest neighbours given a 3D cartesian point (x,y,z)
    NodeList kNearestNeighbours( const Point3& p3, size_t k ) const { return tree_.kNearestNeighbours( p3, k ); }

    /// @brief Find nearest neighbour given a 2D lonlat point (lon,lat)
    Node nearestNeighbour( const Point2& p2 ) const {
        Point3 p3;
        util::Earth::convertSphericalToCartesian( p2, p3 );
        return nearestNeighbour( p3 );
    }

    /// @brief Find nearest neighbour given a 3D cartesian point (x,y,z)
    Node nearestNeighbour( const Point3& p3 ) const { return tree_.nearestNeighbour( p3 ); }

    /// @brief Find all points in within a distance of given radius from a given point (lon,lat)
    NodeList findInSphere( const Point2& p2, double radius ) const {
        Point3 p3;
        util::Earth::convertSphericalToCartesian( p2, p3 );
        return findInSphere( p3, radius );
    }

    /// @brief Find all points in within a distance of given radius from a given point (x,y,z)
    NodeList findInSphere( const Point3& p3, double radius ) const { return tree_.findInSphere( p3, radius ); }

    const eckit_KDTree& tree() const { return tree_; }

private:
    std::vector<typename eckit_KDTree::Value> tmp_;
    eckit_KDTree tree_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
