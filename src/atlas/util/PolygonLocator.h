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

#include <algorithm>
#include <memory>

#include "atlas/library/config.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

///@brief Find polygon that contains a point
///
/// Construction requires a list of polygons.
/// The implementation makes use of a KDTree that holds centroids of all polygons.
/// Once the nearest polygon-centroids (default=4) are found, the corresponding polygons are
/// visited in order of shortest distance, to check if the point is contained within.
class PolygonLocator {
public:
    /// @brief Construct PolygonLocator from shared_ptr of polygons
    PolygonLocator( const std::shared_ptr<const PolygonCoordinates::Vector> polygons,
                    const Projection& projection = Projection() ) :
        shared_polygons_( polygons ), polygons_( *shared_polygons_ ), projection_( projection ) {
        k_ = std::min( k_, polygons_.size() );
        buildKDTree();
    }

    /// @brief Construct PolygonLocator and move polygons inside.
    PolygonLocator( PolygonCoordinates::Vector&& polygons, const Projection& projection = Projection() ) :
        shared_polygons_( std::make_shared<PolygonCoordinates::Vector>( std::move( polygons ) ) ),
        polygons_( *shared_polygons_ ),
        projection_( projection ) {
        k_ = std::min( k_, polygons_.size() );
        buildKDTree();
    }

    /// @brief Construct PolygonLocator using reference to polygons.
    /// !WARNING! polygons should not go out of scope before PolygonLocator
    PolygonLocator( const PolygonCoordinates::Vector& polygons, const Projection& projection = Projection() ) :
        polygons_( polygons ), projection_( projection ) {
        k_ = std::min( k_, polygons_.size() );
        buildKDTree();
    }

    /// @brief find the polygon that holds the point (lon,lat)
    idx_t operator()( const Point2& point ) const {
        const auto found = kdtree_.kNearestNeighbours( point, k_ );
        idx_t partition{-1};
        for ( size_t i = 0; i < found.size(); ++i ) {
            idx_t ii = found[i].payload();
            if ( polygons_[ii].contains( lonlat2xy( point ) ) ) {
                partition = ii;
                break;
            }
        }
        ATLAS_ASSERT( partition >= 0, "Could not find find point in `k` \"nearest\" polygons. Increase `k`?" );
        return partition;
    }

private:
    void buildKDTree() {
        kdtree_.reserve( polygons_.size() );
        for ( idx_t p = 0; p < polygons_.size(); ++p ) {
            kdtree_.insert( xy2lonlat( polygons_[p].centroid() ), p );
        }
        kdtree_.build();
    }

    Point2 lonlat2xy( const Point2& lonlat ) const {
        Point2 xy{lonlat};
        projection_.lonlat2xy( xy.data() );
        return xy;
    }
    Point2 xy2lonlat( const Point2& xy ) const {
        Point2 lonlat{xy};
        projection_.xy2lonlat( lonlat.data() );
        return lonlat;
    }


    std::shared_ptr<const PolygonCoordinates::Vector> shared_polygons_;
    const PolygonCoordinates::Vector& polygons_;
    Projection projection_;
    idx_t k_{4};
    KDTree<idx_t> kdtree_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
