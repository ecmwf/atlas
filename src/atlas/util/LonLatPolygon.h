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

#include <type_traits>

#include "atlas/util/KDTree.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace util {


//------------------------------------------------------------------------------------------------------

/// @brief Implement PolygonCoordinates::contains for a polygon defined in LonLat space.
class LonLatPolygon : public PolygonCoordinates {
private:
    template <typename PointContainer>
    using enable_if_not_polygon = typename std::enable_if<!std::is_base_of<Polygon, PointContainer>::value, int>::type;

public:
    // -- Constructors

    LonLatPolygon( const Polygon&, const atlas::Field& coordinates, bool removeAlignedPoints = true );
    LonLatPolygon( const PartitionPolygon& );

    template <typename PointContainer, enable_if_not_polygon<PointContainer> = 0>
    LonLatPolygon( const PointContainer& points, bool removeAlignedPoints = true );

    /// @brief Point-in-polygon test based on winding number
    /// @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a Polygon</a>
    /// @param[in] P given point
    /// @return if point is in polygon
    bool contains( const Point2& P ) const override;

private:
    PointLonLat centroid_;
    double inner_radius_squared_{0};
    PointLonLat inner_coordinatesMin_;
    PointLonLat inner_coordinatesMax_;
};


//------------------------------------------------------------------------------------------------------

/// @brief Vector of all polygons with functionality to find partition using a KDTree
class LonLatPolygons : public VectorOfAbstract<PolygonCoordinates> {
public:
    LonLatPolygons( const PartitionPolygons& partition_polygons ) {
        reserve( partition_polygons.size() );
        for ( auto& partition_polygon : partition_polygons ) {
            emplace_back( new LonLatPolygon( partition_polygon ) );
        }
        search_.reserve( size() );
        for ( idx_t p = 0; p < size(); ++p ) {
            search_.insert( this->at( p ).centroid(), p );
        }
        search_.build();
        k_ = std::min( k_, size() );
    }

    /// @brief find the partition that holds the point (lon,lat)
    idx_t findPartition( const Point2& point ) {
        const auto found = search_.kNearestNeighbours( point, k_ );
        idx_t partition{-1};
        for ( size_t i = 0; i < found.size(); ++i ) {
            idx_t ii = found[i].payload();
            if ( this->at( ii ).contains( point ) ) {
                partition = ii;
                break;
            }
        }
        ASSERT( partition >= 0 );
        return partition;
    }

private:
    KDTree<idx_t> search_;
    idx_t k_{4};
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
