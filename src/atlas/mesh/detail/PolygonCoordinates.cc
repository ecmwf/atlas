/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/PolygonCoordinates.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/detail/Polygon.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace mesh {
namespace detail {


//------------------------------------------------------------------------------------------------------


namespace {


double cross_product_analog(const PointLonLat& A, const PointLonLat& B, const PointLonLat& C) {
  return (A[LON] - C[LON]) * (B[LAT] - C[LAT])
       - (A[LAT] - C[LAT]) * (B[LON] - C[LON]);
}


}  // (anonymous)


//------------------------------------------------------------------------------------------------------


PolygonCoordinates::PolygonCoordinates(
        const Polygon& poly,
        const atlas::Field& lonlat,
        bool includesNorthPole,
        bool includesSouthPole,
        bool removeAlignedPoints ) :
    includesNorthPole_(includesNorthPole),
    includesSouthPole_(includesSouthPole) {

    ASSERT(poly.size() > 2);
    ASSERT(poly.front() == poly.back());

    // Point coordinates
    // - use a bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included

    coordinates_.clear();
    coordinates_.reserve(poly.size());

    auto ll = array::make_view< double, 2 >(lonlat);
    coordinatesMin_ = PointLonLat(ll(poly[0], LON), ll(poly[0], LAT));
    coordinatesMax_ = coordinatesMin_;

    for (size_t i = 0; i < poly.size(); ++i) {

        PointLonLat A(ll(poly[i], LON), ll(poly[i], LAT));
        coordinatesMin_ = PointLonLat::componentsMin(coordinatesMin_, A);
        coordinatesMax_ = PointLonLat::componentsMax(coordinatesMax_, A);

        // if new point is aligned with existing edge (cross product ~= 0) make the edge longer
        if ((coordinates_.size() >= 2) && removeAlignedPoints) {
            const PointLonLat& B = coordinates_.back();
            const PointLonLat& C = coordinates_[coordinates_.size() - 2];
            if (eckit::types::is_approximately_equal( 0., cross_product_analog(A, B, C) )) {
                coordinates_.back() = A;
                continue;
            }
        }

        coordinates_.push_back(A);
    }

    ASSERT(coordinates_.size() == poly.size());
}


PolygonCoordinates::PolygonCoordinates(
        const std::vector<PointLonLat>& points,
        bool includesNorthPole,
        bool includesSouthPole ) :
    coordinates_(points),
    includesNorthPole_(includesNorthPole),
    includesSouthPole_(includesSouthPole) {

    ASSERT(coordinates_.size() > 2);
    ASSERT(eckit::geometry::points_equal(coordinates_.front(), coordinates_.back()));

    coordinatesMin_ = coordinates_.front();
    coordinatesMax_ = coordinatesMin_;
    for (const PointLonLat& P : coordinates_) {
        coordinatesMin_ = PointLonLat::componentsMin(coordinatesMin_, P);
        coordinatesMax_ = PointLonLat::componentsMax(coordinatesMax_, P);
    }
}


bool PolygonCoordinates::containsPointInLonLatGeometry(const PointLonLat& P) const {
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

            // check point-edge side and direction, using 2D-analog cross-product;
            // tests if P is left|on|right of a directed A-B infinite line, by intersecting either:
            // - "up" on upward crossing & P left of edge, or
            // - "down" on downward crossing & P right of edge
            const bool APB = (A[LAT] <= P[LAT] && P[LAT] < B[LAT]);
            const bool BPA = (B[LAT] <= P[LAT] && P[LAT] < A[LAT]);

            if (APB != BPA) {
                const double side = cross_product_analog(P, A, B);
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


bool PolygonCoordinates::containsPointInSphericalGeometry(const PointLonLat& P) const {
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

            // test if P is on/above/below of a great circle containing A,B
            const bool APB = (A[LON] <= P[LON] && P[LON] < B[LON]);
            const bool BPA = (B[LON] <= P[LON] && P[LON] < A[LON]);

            if (APB != BPA) {
                PointLonLat p(P[LON], std::numeric_limits<double>::quiet_NaN());
                util::Earth::greatCircleLatitudeGivenLongitude(A, B, p);

                ASSERT(!std::isnan(p.lat()));
                if (eckit::types::is_approximately_equal(P[LAT], p.lat())) {
                    return true;
                }

                wn += (P[LAT] > p.lat() ? -1:1) * (APB ? -1:1);
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

