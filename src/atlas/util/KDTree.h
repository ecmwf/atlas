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

#include "eckit/config/Configuration.h"

#include "atlas/util/Geometry.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/detail/KDTree.h"

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
///     auto neighbours = search.closestPoints( PointLonLat{180., 45.}, k ).payloads();
/// @endcode
/// The variable `neighbours` is now a container of indices (the payloads) of the 4 nearest points

template <typename PayloadT, typename PointT = Point3>
class KDTree : public ObjectHandle<detail::KDTreeBase<PayloadT, PointT>> {
public:
    using Handle         = ObjectHandle<detail::KDTreeBase<PayloadT, PointT>>;
    using Implementation = typename Handle::Implementation;
    using Point          = typename Implementation::Point;
    using Payload        = typename Implementation::Payload;
    using PayloadList    = typename Implementation::PayloadList;
    using Value          = typename Implementation::Value;
    using ValueList      = typename Implementation::ValueList;

    using Handle::get;
    using Handle::Handle;

public:
    //--------------------------------------------------------------------------------------
    // Constructors

    /// @brief Construct an empty kd-tree with default geometry (Earth)
    KDTree(): Handle(new detail::KDTreeMemory<Payload, Point>()) {}

    /// @brief Construct an empty kd-tree with custom geometry
    KDTree(const Geometry& geometry): Handle(new detail::KDTreeMemory<Payload, Point>(geometry)) {}

    /// @brief Construct an empty kd-tree with custom geometry
    KDTree(const eckit::Configuration& config):
        KDTree(Geometry(config.getString("geometry","Earth"))) {}

    /// @brief Construct a shared kd-tree with default geometry (Earth)
    template <typename Tree>
    KDTree(const std::shared_ptr<Tree>& kdtree): Handle(new detail::KDTree_eckit<Tree, Payload, Point>(kdtree)) {}

    /// @brief Construct a shared kd-tree with custom geometry
    template <typename Tree>
    KDTree(const std::shared_ptr<Tree>& kdtree, const Geometry& geometry):
        Handle(new detail::KDTree_eckit<Tree, Payload, Point>(kdtree, geometry)) {}

    //--------------------------------------------------------------------------------------
    // Methods to build the KDTree

    /// @brief Reserve memory for building the kd-tree in one shot (optional, at cost of extra memory)
    void reserve(idx_t size) { get()->reserve(size); }

    /// @brief Insert spherical point (lon,lat) or 3D cartesian point (x,y,z)
    /// @warning If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    template <typename Point>
    void insert(const Point& p, const Payload& payload) {
        get()->insert(p, payload);
    }

    /// @brief Insert kd-tree value in tree
    /// @warning If memory has been reserved with reserve(), insertion will be delayed until build() is called.
    void insert(const Value& value) { get()->insert(value); }

    /// @brief Build the kd-tree in one shot, if memory has been reserved.
    /// This will need to be called before all search functions like closestPoints().
    /// @post The KDTree is ready to be used
    void build() { get()->build(); }

    /// @brief Build with spherical points (lon,lat) where longitudes, latitudes, and payloads are separate containers.
    /// Memory will be reserved with reserve() to match the size
    /// @post The KDTree is ready to be used
    template <typename Longitudes, typename Latitudes, typename Payloads>
    void build(const Longitudes& longitudes, const Latitudes& latitudes, const Payloads& payloads) {
        get()->build(longitudes, latitudes, payloads);
    }

    /// @brief Build with spherical points (lon,lat) given separate iterator ranges for longitudes, latitudes, and payloads.
    /// Memory will be reserved with reserve() to match the size
    /// @post The KDTree is ready to be used
    template <typename LongitudesIterator, typename LatitudesIterator, typename PayloadsIterator>
    void build(const LongitudesIterator& longitudes_begin, const LongitudesIterator& longitudes_end,
               const LatitudesIterator& latitudes_begin, const LatitudesIterator& latitudes_end,
               const PayloadsIterator& payloads_begin, const PayloadsIterator& payloads_end) {
        get()->build(longitudes_begin, longitudes_end, latitudes_begin, latitudes_end, payloads_begin, payloads_end);
    }

    /// @brief Build with templated points. Points can be either 2D (lon,lat) or 3D (x,y,z)
    /// Memory will be reserved with reserve() to match the size
    /// @post The KDTree is ready to be used
    template <typename Points, typename Payloads>
    void build(const Points& points, const Payloads& payloads) {
        get()->build(points, payloads);
    }

    /// @brief Build with spherical points (lon,lat) given separate iterator ranges for longitudes, latitudes, and payloads.
    /// Memory will be reserved with reserve() to match the size
    /// @post The KDTree is ready to be used
    template <typename PointIterator, typename PayloadsIterator>
    void build(const PointIterator& points_begin, const PointIterator& points_end,
               const PayloadsIterator& payloads_begin, const PayloadsIterator& payloads_end) {
        get()->build(points_begin, points_end, payloads_begin, payloads_end);
    }

    /// @brief Build with vector of Value
    /// @post The KDTree is ready to be used
    void build(const std::vector<Value>& values) { get()->build(values); }

    //--------------------------------------------------------------------------------------
    // Methods to access the KDTree

    bool empty() const { return get()->empty(); }

    size_t size() const { return get()->size(); }

    size_t footprint() const { return get()->footprint(); }

    /// @brief Find k closest points given a 3D cartesian point (x,y,z) or 2D lonlat point(lon,lat)
    template <typename Point>
    ValueList closestPoints(const Point& p, size_t k) const {
        return get()->closestPoints(p, k);
    }

    /// @brief Find closest point given a 3D cartesian point (x,y,z) or 2D lonlat point(lon,lat)
    template <typename Point>
    Value closestPoint(const Point& p) const {
        return get()->closestPoint(p);
    }

    /// @brief Find all points within a distance of given radius from a given point 3D cartesian point (x,y,z)
    /// or a 2D (lon,lat) point
    template <typename Point>
    ValueList closestPointsWithinRadius(const Point& p, double radius) const {
        return get()->closestPointsWithinRadius(p, radius);
    }

    /// @brief Return geometry used to convert (lon,lat) to (x,y,z) coordinates
    const Geometry& geometry() const { return get()->geometry(); }
};

//------------------------------------------------------------------------------------------------------

using IndexKDTree2D = KDTree<idx_t, Point2>;  // 2D search: implementation is using 2D points only
using IndexKDTree3D = KDTree<idx_t, Point3>;  // 3D search: lonlat (2D) to xyz (3D) conversion is done internally
using IndexKDTree   = IndexKDTree3D;

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
IndexKDTree::Implementation* atlas__IndexKDTree__new();
IndexKDTree::Implementation* atlas__IndexKDTree__new_geometry(const Geometry::Implementation* geometry);
void atlas__IndexKDTree__delete(IndexKDTree::Implementation* This);
void atlas__IndexKDTree__reserve(IndexKDTree::Implementation* This, const idx_t size);
void atlas__IndexKDTree__insert(IndexKDTree::Implementation* This, const double lon, const double lat,
                                const idx_t index);
void atlas__IndexKDTree__build(IndexKDTree::Implementation* This);
void atlas__IndexKDTree__closestPoints(const IndexKDTree::Implementation* This, const double plon, const double plat,
                                       const size_t k, double*& lon, double*& lat, idx_t*& indices, double*& distances);
void atlas__IndexKDTree__closestPoint(const IndexKDTree::Implementation* This, const double plon, const double plat,
                                      double& lon, double& lat, idx_t& index, double& distance);
void atlas__IndexKDTree__closestPointsWithinRadius(const IndexKDTree::Implementation* This, const double plon,
                                                   const double plat, const double radius, size_t& k, double*& lon,
                                                   double*& lat, idx_t*& indices, double*& distances);
const Geometry::Implementation* atlas__IndexKDTree__geometry(const IndexKDTree::Implementation* This);
int atlas__IndexKDTree__empty(const IndexKDTree::Implementation* This);
idx_t atlas__IndexKDTree__size(const IndexKDTree::Implementation* This);
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
