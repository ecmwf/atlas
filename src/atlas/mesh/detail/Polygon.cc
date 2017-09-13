/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/Polygon.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {
namespace detail {

//------------------------------------------------------------------------------------------------------


namespace {


double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Ay - Cy) * (Bx - Cx);
}


/*
 * Tests if a given point is left|on|right of an infinite line.
 * @input P point to test
 * @input A, B points on infinite line
 * @return >0/=0/<0 for P left|on|right of the infinite line
 */
double point_on_which_side(const PointLonLat& P, const PointLonLat& A, const PointLonLat& B) {
    return dot_sign( P[LON], P[LAT],
                     A[LON], A[LAT],
                     B[LON], B[LAT] );
}


}  // (anonymous)


//------------------------------------------------------------------------------------------------------


Polygon::Polygon() {
}


Polygon::Polygon(const Polygon::edge_set_t& edges) {

    // get external edges by attempting to remove reversed edges, if any
    edge_set_t extEdges;
    for (const edge_t& e : edges) {
        if (!extEdges.erase(e.reverse())) {
            extEdges.insert(e);
        }
    }
    ASSERT(extEdges.size() >= 2);

    // set one polygon cycle, by picking next edge with first node same as second node of last edge
    clear();
    reserve(extEdges.size() + 1);

    push_back(extEdges.begin()->first);
    for (edge_set_t::iterator e = extEdges.begin(); e != extEdges.end() && e->first == back();
         e = extEdges.lower_bound(edge_t(back(), std::numeric_limits< idx_t >::min()))) {
        push_back(e->second);
        extEdges.erase(*e);
    }
    ASSERT(front() == back());

    // exhaust remaining edges, should represent additional cycles, if any
    while (!extEdges.empty()) {
        operator+=(Polygon(extEdges));
    }
}


Polygon::operator bool() const {
    return !Polygon::empty();
}


Polygon& Polygon::operator+=(const Polygon& other) {

    if (empty()) {
        return operator=(other);
    }

    // polygon can have multiple cycles, but must be connected graphs
    // Note: a 'cycle' is handled by repeating the indices, excluding (repeated) last index
    const difference_type N = difference_type(other.size()) - 1;

    container_t cycle;
    cycle.reserve(2 * size_t(N));
    cycle.insert(cycle.end(), other.begin(), other.end() - 1);
    cycle.insert(cycle.end(), other.begin(), other.end() - 1);

    for (const_iterator c = cycle.begin(); c != cycle.begin() + N; ++c) {
        iterator here = std::find(begin(), end(), *c);
        if (here != end()) {
            insert(here, c, c + N);
            return *this;
        }
    }

    throw eckit::AssertionFailed("Polygon: could not merge polygons, they are not connected", Here());
}


void Polygon::print(std::ostream& s) const {
    char z = '{';
    for (auto n : static_cast<const container_t&>(*this)) {
        s << z << n;
        z = ',';
    }
    s << '}';
}


bool Polygon::containsPointInLonLatGeometry(const PointLonLat& P) const {
    ASSERT(coordinates_.size() >= 2);

    // check first bounding box
    if (coordinatesMin_[LON] <= P[LON] && P[LON] < coordinatesMax_[LON]
     && coordinatesMin_[LAT] <= P[LAT] && P[LAT] < coordinatesMax_[LAT]) {

        // winding number
        int wn = 0;

        // loop on polygon edges
        for (size_t i = 1; i < coordinates_.size(); ++i) {
            const PointLonLat& A = coordinates_[i-1];
            const PointLonLat& B = coordinates_[ i ];

            // intersect either:
            // - "up" on upward crossing & P left of edge, or
            // - "down" on downward crossing & P right of edge
            const bool APB = (A[LAT] <= P[LAT] && P[LAT] < B[LAT]);
            const bool BPA = (B[LAT] <= P[LAT] && P[LAT] < A[LAT]);

            if (APB || BPA) {
                const double side = point_on_which_side(P, A, B);
                if (APB && side > 0) {
                    ++wn;
                } else if (BPA && side < 0) {
                    --wn;
                }
            }
        }

        // wn == 0 only when P is outside
        return wn != 0;
    }

    return ((includesNorthPole_ && P[LAT] >= coordinatesMax_[LAT])
         || (includesSouthPole_ && P[LAT] <  coordinatesMin_[LAT]));
}


bool Polygon::containsPointInSphericalGeometry(const PointLonLat& P) const {
    ASSERT(coordinates_.size() >= 2);

    // check first bounding box
    if (coordinatesMin_[LON] <= P[LON] && P[LON] < coordinatesMax_[LON]
     && coordinatesMin_[LAT] <= P[LAT] && P[LAT] < coordinatesMax_[LAT]) {

        // winding number
        int wn = 0;

        // loop on polygon edges
        for (size_t i = 1; i < coordinates_.size(); ++i) {
            const PointLonLat& A = coordinates_[i-1];
            const PointLonLat& B = coordinates_[ i ];

            // intersect either:
            // - "left" on left-wards crossing & P above edge, or
            // - "right" on right-wards crossing & P below edge
            const bool APB = (A[LON] <= P[LON] && P[LON] < B[LON]);
            const bool BPA = (B[LON] <= P[LON] && P[LON] < A[LON]);

            if (APB || BPA) {
                const double side = point_on_which_side(P, A, B);
                if (eckit::types::is_approximately_equal(side, 0.)) {
                    return true;
                } else if (APB && side > 0) {
                    ++wn;
                } else if (BPA && side < 0) {
                    --wn;
                }
            }
        }

        // wn == 0 only when P is outside
        return wn != 0;
    }

    return ((includesNorthPole_ && P[LAT] >= coordinatesMax_[LAT])
         || (includesSouthPole_ && P[LAT] <  coordinatesMin_[LAT]));
}


//------------------------------------------------------------------------------------------------------

}  // detail
}  // mesh
}  // atlas

