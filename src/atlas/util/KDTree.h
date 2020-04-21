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
#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

/// @brief k-dimensional tree constructable both with 2D (lon,lat) points as with 3D (x,y,z) points
///
/// The implementation is based on eckit::KDTreeMemory with 3D (x,y,z) points.
/// 2D points (lon,lat) are converted when needed to 3D during insertion, and during search, so that
/// a search always happens with 3D cartesian points.
///
/// ### Example:
///
/// Construct a KDTree, given a `list_of_lonlat_points` (e.g. `std::vector<PointLonLat>` or other container)
/// A payload given is the index in this list.
/// @code{.cpp}
///     KDTree<idx_t> search;
///     search.reserve( list_of_lonlat_points.size() );
///     idx_t n{0};
///     for ( auto& point : list_of_lonlat_points ) {
///         search.insert( point, n++ );
///     }
///     search.build();
/// @endcode
/// We can now do e.g. a search for the nearest 4 neighbours (k=4) sorted by shortest distance
/// @code{.cpp}
///     idx_t k = 4;
///     auto neighbours = search.kNearestNeighbours( PointLonLat{180., 45.}, k ).payloads();
/// @endcode
/// The variable `neighbours` is now a container of indices (the payloads) of the 4 nearest points
template <typename PayloadT>
class KDTree {
private:
    struct eckit_KDTreeTraits {
        using Point   = Point3;
        using Payload = PayloadT;
    };
    using eckit_KDTree = eckit::KDTreeMemory<eckit_KDTreeTraits>;

public:
    using Payload     = PayloadT;
    using Node        = typename eckit_KDTree::NodeInfo;
    using PayLoadList = std::vector<Payload>;

    class NodeList : public eckit_KDTree::NodeList {
        using Base = typename eckit_KDTree::NodeList;

    public:
        NodeList( typename eckit_KDTree::NodeList&& list ) : Base( std::move( list ) ) {}
        PayLoadList payloads() const {
            PayLoadList list;
            list.reserve( Base::size() );
            for ( auto& item : *this ) {
                list.emplace_back( item.payload() );
            }
            return list;
        }
    };

    /// @brief Reserve memory for building the kdtree in one shot (optional, at cost of extra memory)
    void reserve( idx_t size );

    /// @brief Insert spherical point (lon,lat)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void insert( const Point2& p2, const Payload& payload );

    /// @brief Insert 3D cartesian point (x,y,z)
    /// If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void insert( const Point3& p3, const Payload& payload );

    /// @brief Build the kd-tree in one shot, if memory has been reserved.
    /// This will need to be called before all search functions like kNearestNeighbours().
    void build();

    /// @brief Find k nearest neighbour given a 2D lonlat point (lon,lat)
    NodeList kNearestNeighbours( const Point2& p2, size_t k ) const;

    /// @brief Find k nearest neighbours given a 3D cartesian point (x,y,z)
    NodeList kNearestNeighbours( const Point3& p3, size_t k ) const;

    /// @brief Find nearest neighbour given a 2D lonlat point (lon,lat)
    Node nearestNeighbour( const Point2& p2 ) const;

    /// @brief Find nearest neighbour given a 3D cartesian point (x,y,z)
    Node nearestNeighbour( const Point3& p3 ) const;

    /// @brief Find all points in within a distance of given radius from a given point (lon,lat)
    NodeList findInSphere( const Point2& p2, double radius ) const;

    /// @brief Find all points in within a distance of given radius from a given point (x,y,z)
    NodeList findInSphere( const Point3& p3, double radius ) const;

    const eckit_KDTree& tree() const { return tree_; }

private:
    void assert_built() const;

private:
    std::vector<typename eckit_KDTree::Value> tmp_;
    mutable eckit_KDTree tree_;  // mutable because its member functions are non-const...
};

//------------------------------------------------------------------------------------------------------

template <typename PayloadT>
void KDTree<PayloadT>::reserve( idx_t size ) {
    tmp_.reserve( size );
}

template <typename PayloadT>
void KDTree<PayloadT>::insert( const Point2& p2, const Payload& payload ) {
    Point3 p3;
    util::Earth::convertSphericalToCartesian( p2, p3 );
    insert( p3, payload );
}

template <typename PayloadT>
void KDTree<PayloadT>::insert( const Point3& p3, const Payload& payload ) {
    if ( tmp_.capacity() ) {
        tmp_.emplace_back( p3, payload );
    }
    else {
        tree_.insert( {p3, payload} );
    }
}

template <typename PayloadT>
void KDTree<PayloadT>::build() {
    if ( tmp_.size() ) {
        tree_.build( tmp_.begin(), tmp_.end() );
        tmp_.clear();
        tmp_.shrink_to_fit();
    }
}

template <typename PayloadT>
typename KDTree<PayloadT>::NodeList KDTree<PayloadT>::kNearestNeighbours( const Point2& p2, size_t k ) const {
    Point3 p3;
    util::Earth::convertSphericalToCartesian( p2, p3 );
    return kNearestNeighbours( p3, k );
}

template <typename PayloadT>
typename KDTree<PayloadT>::NodeList KDTree<PayloadT>::kNearestNeighbours( const Point3& p3, size_t k ) const {
    assert_built();
    return tree_.kNearestNeighbours( p3, k );
}

template <typename PayloadT>
typename KDTree<PayloadT>::Node KDTree<PayloadT>::nearestNeighbour( const Point2& p2 ) const {
    Point3 p3;
    util::Earth::convertSphericalToCartesian( p2, p3 );
    return nearestNeighbour( p3 );
}

template <typename PayloadT>
typename KDTree<PayloadT>::Node KDTree<PayloadT>::nearestNeighbour( const Point3& p3 ) const {
    assert_built();
    return tree_.nearestNeighbour( p3 );
}

template <typename PayloadT>
typename KDTree<PayloadT>::NodeList KDTree<PayloadT>::findInSphere( const Point2& p2, double radius ) const {
    Point3 p3;
    util::Earth::convertSphericalToCartesian( p2, p3 );
    return findInSphere( p3, radius );
}

template <typename PayloadT>
typename KDTree<PayloadT>::NodeList KDTree<PayloadT>::findInSphere( const Point3& p3, double radius ) const {
    assert_built();
    return tree_.findInSphere( p3, radius );
}

template <typename PayloadT>
void KDTree<PayloadT>::assert_built() const {
    if ( tmp_.capacity() ) {
        throw_AssertionFailed( "KDTree was used before calling build()" );
    }
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
