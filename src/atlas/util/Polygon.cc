/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace util {

namespace {

//------------------------------------------------------------------------------------------------------

double cross_product_analog( const PointLonLat& A, const PointLonLat& B, const PointLonLat& C ) {
    return ( A.lon() - C.lon() ) * ( B.lat() - C.lat() ) - ( A.lat() - C.lat() ) * ( B.lon() - C.lon() );
}

}  // namespace

//------------------------------------------------------------------------------------------------------

Polygon::Polygon() {}

Polygon::Polygon( const Polygon::edge_set_t& edges ) {
    ATLAS_TRACE();
    // get external edges by attempting to remove reversed edges, if any
    edge_set_t extEdges;
    for ( const edge_t& e : edges ) {
        if ( !extEdges.erase( e.reverse() ) ) {
            extEdges.insert( e );
        }
    }
    ATLAS_ASSERT( extEdges.size() >= 2 );

    // set one polygon cycle, by picking next edge with first node same as second
    // node of last edge
    clear();
    reserve( extEdges.size() + 1 );

    push_back( extEdges.begin()->first );
    for ( edge_set_t::iterator e = extEdges.begin(); e != extEdges.end() && e->first == back();
          e                      = extEdges.lower_bound( edge_t( back(), std::numeric_limits<idx_t>::min() ) ) ) {
        push_back( e->second );
        extEdges.erase( *e );
    }
    ATLAS_ASSERT( front() == back() );

    // exhaust remaining edges, should represent additional cycles, if any
    while ( !extEdges.empty() ) {
        operator+=( Polygon( extEdges ) );
    }
}

Polygon::operator bool() const {
    return !Polygon::empty();
}

Polygon& Polygon::operator+=( const Polygon& other ) {
    if ( empty() ) {
        return operator=( other );
    }

    // polygon can have multiple cycles, but must be connected graphs
    // Note: a 'cycle' is handled by repeating the indices, excluding (repeated)
    // last index
    ATLAS_ASSERT( other.front() == other.back() );
    const difference_type N = difference_type( other.size() ) - 1;

    container_t cycle( 2 * size_t( N ) );
    std::copy( other.begin(), other.begin() + N, cycle.begin() );
    std::copy( other.begin(), other.begin() + N, cycle.begin() + N );

    for ( const_iterator c = cycle.begin(); c != cycle.begin() + N; ++c ) {
        iterator here = std::find( begin(), end(), *c );
        if ( here != end() ) {
            insert( here, c, c + N );
            return *this;
        }
    }

    throw_AssertionFailed( "Polygon: could not merge polygons, they are not connected", Here() );
}

void Polygon::print( std::ostream& s ) const {
    char z = '{';
    for ( auto n : static_cast<const container_t&>( *this ) ) {
        s << z << n;
        z = ',';
    }
    s << '}';
}

//------------------------------------------------------------------------------------------------------

PolygonCoordinates::PolygonCoordinates( const Polygon& poly, const atlas::Field& lonlat, bool removeAlignedPoints ) {
    ATLAS_ASSERT( poly.size() > 2 );
    ATLAS_ASSERT( poly.front() == poly.back() );

    // Point coordinates
    // - use a bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included

    coordinates_.clear();
    coordinates_.reserve( poly.size() );

    auto ll         = array::make_view<double, 2>( lonlat );
    coordinatesMin_ = PointLonLat( ll( poly[0], LON ), ll( poly[0], LAT ) );
    coordinatesMax_ = coordinatesMin_;

    size_t nb_removed_points_due_to_alignment = 0;

    for ( size_t i = 0; i < poly.size(); ++i ) {
        PointLonLat A( ll( poly[i], LON ), ll( poly[i], LAT ) );
        coordinatesMin_ = PointLonLat::componentsMin( coordinatesMin_, A );
        coordinatesMax_ = PointLonLat::componentsMax( coordinatesMax_, A );

        // if new point is aligned with existing edge (cross product ~= 0) make the
        // edge longer
        if ( ( coordinates_.size() >= 2 ) && removeAlignedPoints ) {
            const PointLonLat& B = coordinates_.back();
            const PointLonLat& C = coordinates_[coordinates_.size() - 2];
            if ( eckit::types::is_approximately_equal( 0., cross_product_analog( A, B, C ) ) ) {
                coordinates_.back() = A;
                ++nb_removed_points_due_to_alignment;
                continue;
            }
        }

        coordinates_.push_back( A );
    }

    ATLAS_ASSERT( coordinates_.size() == poly.size() - nb_removed_points_due_to_alignment );
}

PolygonCoordinates::PolygonCoordinates( const std::vector<PointLonLat>& points ) : coordinates_( points ) {
    ATLAS_ASSERT( coordinates_.size() > 2 );
    ATLAS_ASSERT( eckit::geometry::points_equal( coordinates_.front(), coordinates_.back() ) );

    coordinatesMin_ = coordinates_.front();
    coordinatesMax_ = coordinatesMin_;
    for ( const PointLonLat& P : coordinates_ ) {
        coordinatesMin_ = PointLonLat::componentsMin( coordinatesMin_, P );
        coordinatesMax_ = PointLonLat::componentsMax( coordinatesMax_, P );
    }
}

PolygonCoordinates::~PolygonCoordinates() {}

const PointLonLat& PolygonCoordinates::coordinatesMax() const {
    return coordinatesMax_;
}

const PointLonLat& PolygonCoordinates::coordinatesMin() const {
    return coordinatesMin_;
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
