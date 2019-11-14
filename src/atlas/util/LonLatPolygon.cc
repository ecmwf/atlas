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

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatPolygon.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

namespace {

double cross_product_analog( const PointLonLat& A, const PointLonLat& B, const PointLonLat& C ) {
    return ( A.lon() - C.lon() ) * ( B.lat() - C.lat() ) - ( A.lat() - C.lat() ) * ( B.lon() - C.lon() );
}


PointLonLat compute_centroid(const std::vector<PointLonLat>& points) {

    ATLAS_ASSERT( eckit::geometry::points_equal( points.front(), points.back() ) );

    PointLonLat centroid = {0, 0};
    double signed_area = 0.;
    double a = 0.;  // Partial signed area

    for (size_t i=0; i<points.size()-1; ++i) {
        const PointLonLat& p0 = points[i];
        const PointLonLat& p1 = points[i+1];
        a = p0[0]*p1[1] - p1[0]*p0[1];
        signed_area += a;
        centroid[0] += (p0[0] + p1[0])*a;
        centroid[1] += (p0[1] + p1[1])*a;
    }
    signed_area *= 0.5;
    centroid[0] /= (6.*signed_area);
    centroid[1] /= (6.*signed_area);

    return centroid;
}

double compute_inner_radius_squared(const std::vector<PointLonLat>& points, const PointLonLat& centroid) {
  auto distance2 = [](const PointLonLat& a, const PointLonLat& b) {
    double dx = (a[0]-b[0]);
    double dy = (a[1]-b[1]);
    return dx*dx + dy*dy;
  };
  auto dot = [](const PointLonLat& a, const PointLonLat& b) {
    return a[0]*b[0]+a[1]*b[1];
  };
  double R2 = std::numeric_limits<double>::max();
  PointLonLat projection;
  for( size_t i=0; i<points.size()-1; ++i ) {
    double d2 = distance2(points[i],points[i+1]);
    double t = std::max(0.,std::min(1.,dot(centroid-points[i],points[i+1]-points[i])/d2));
    projection[0] = points[i][0] + t * (points[i+1][0]-points[i][0]);
    projection[1] = points[i][1] + t * (points[i+1][1]-points[i][1]);
    R2 = std::min( R2, distance2(projection,centroid) );
    //Log::info() << "Segment " << points[i] << " - " << points[i+1] << " :  projection = " << projection << "   \t distance = "  << std::sqrt(distance2(projection,centroid) ) << std::endl;
  }
  return R2;
}

}  // namespace

//------------------------------------------------------------------------------------------------------

LonLatPolygon::LonLatPolygon( const Polygon& poly, const atlas::Field& lonlat, bool removeAlignedPoints ) :
    PolygonCoordinates( poly, lonlat, removeAlignedPoints ) {

    centroid_ = compute_centroid( coordinates_ );
    inner_radius_squared_ = compute_inner_radius_squared( coordinates_, centroid_ );

    }

LonLatPolygon::LonLatPolygon( const std::vector<PointLonLat>& points, bool removeAlignedPoints ) :
    PolygonCoordinates( points, removeAlignedPoints ) {

    centroid_ = compute_centroid( coordinates_ );
    inner_radius_squared_ = compute_inner_radius_squared( coordinates_, centroid_ );

    ATLAS_ASSERT( contains(centroid_) );
}

LonLatPolygon::LonLatPolygon( const std::vector<PointLonLat>& points, const std::vector<PointLonLat>& inner_bounding_box, bool removeAlignedPoints ) : 
    PolygonCoordinates( points, removeAlignedPoints ) {

    inner_coordinatesMin_ = PointLonLat( std::numeric_limits<double>::max(), std::numeric_limits<double>::max() );
    inner_coordinatesMax_ = PointLonLat( std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() );
    for ( size_t i = 0; i < inner_bounding_box.size(); ++i ) {
        inner_coordinatesMin_ = PointLonLat::componentsMin( inner_coordinatesMin_, inner_bounding_box[i] );
        inner_coordinatesMax_ = PointLonLat::componentsMax( inner_coordinatesMax_, inner_bounding_box[i] );
    }
}

bool LonLatPolygon::contains( const PointLonLat& P ) const {

    // check first bounding box
    if ( coordinatesMax_.lat() < P.lat() ||
         P.lat() < coordinatesMin_.lat() || coordinatesMax_.lon() < P.lon() || P.lon() < coordinatesMin_.lon() ) {
        return false;
    }

    auto distance2 = [](const PointLonLat& p, const PointLonLat& centroid) {
      double dx = (p[0]-centroid[0]);
      double dy = (p[1]-centroid[1]);
      return dx*dx + dy*dy;
    };

    if( inner_radius_squared_ == 0 ) { // check inner bounding box
        if ( inner_coordinatesMin_.lon() <= P.lon() &&
             inner_coordinatesMax_.lon() >= P.lon() &&
             inner_coordinatesMin_.lat() <= P.lat() &&
             inner_coordinatesMax_.lat() >= P.lat() ) {
             return true;
        }
    }
    else {     // check inscribed circle
        if( distance2(P,centroid_) < inner_radius_squared_ ) {
            return true;
        }
    }

    // winding number
    int wn = 0;

    // loop on polygon edges
    for ( size_t i = 1; i < coordinates_.size(); ++i ) {
        const PointLonLat& A = coordinates_[i - 1];
        const PointLonLat& B = coordinates_[i];

        // check point-edge side and direction, using 2D-analog cross-product;
        // tests if P is left|on|right of a directed A-B infinite line, by
        // intersecting either:
        // - "up" on upward crossing & P left of edge, or
        // - "down" on downward crossing & P right of edge
        const bool APB = ( A.lat() <= P.lat() && P.lat() <= B.lat() );
        const bool BPA = ( B.lat() <= P.lat() && P.lat() <= A.lat() );

        if ( APB != BPA ) {
            const double side = cross_product_analog( P, A, B );
            if ( APB && side > 0 ) {
                ++wn;
            }
            else if ( BPA && side < 0 ) {
                --wn;
            }
        }
    }

    // wn == 0 only when P is outside
    return wn != 0;
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
