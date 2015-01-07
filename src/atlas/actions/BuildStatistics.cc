/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>
#include <cmath>
#include <set>
#include <limits>
#include <iostream>
#include <algorithm>    // std::sort

#include "eckit/geometry/Point2.h"
#include "atlas/atlas.h"
#include "atlas/mpl/Checksum.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/Parameters.h"
#include "atlas/Util.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"

using eckit::geometry::LLPoint2;
namespace atlas {
namespace actions {

namespace {

/// @brief The usual PI/180 constant
static const double DEG_TO_RAD = M_PI/180.;

/** @brief Computes the arc, in radian, between two positions.
  *
  * The result is equal to <code>Distance(from,to)/EARTH_RADIUS_IN_METERS</code>
  *    <code>= 2*asin(sqrt(h(d/EARTH_RADIUS_IN_METERS )))</code>
  *
  * where:<ul>
  *    <li>d is the distance in meters between 'from' and 'to' positions.</li>
  *    <li>h is the haversine function: <code>h(x)=sin²(x/2)</code></li>
  * </ul>
  *
  * The haversine formula gives:
  *    <code>h(d/R) = h(from.lat-to.lat)+h(from.lon-to.lon)+cos(from.lat)*cos(to.lat)</code>
  *
  * @sa http://en.wikipedia.org/wiki/Law_of_haversines
  */
double arc_in_rad(const LLPoint2& from, const LLPoint2& to) {
    double lat_arc = (from.lat() - to.lat()) * DEG_TO_RAD;
    double lon_arc = (from.lon() - to.lon()) * DEG_TO_RAD;
    double lat_H = std::sin(lat_arc * 0.5);
    lat_H *= lat_H;
    double lon_H = std::sin(lon_arc * 0.5);
    lon_H *= lon_H;
    double tmp = std::cos(from.lat()*DEG_TO_RAD) * std::cos(to.lat()*DEG_TO_RAD);
    return 2.0 * std::asin(std::sqrt(lat_H + tmp*lon_H));
}

}

void build_statistics( Mesh& mesh )
{
  FunctionSpace& nodes = mesh.function_space( "nodes" );
  ArrayView<double,2> coords ( nodes.field( "coordinates"    ) );
  ArrayView<double,1> dist_lon_deg ( nodes.create_field<double>("dist_lon_deg",1) );

  if( mesh.has_function_space("edges") )
  {
    FunctionSpace& edges = mesh.function_space( "edges" );
    if( !edges.has_field("length_deg") )
      edges.create_field<double>("length_deg",1);
    ArrayView<double,1> dist ( edges.field("length_deg") );

    IndexView<int,2> edge_nodes ( edges.field("nodes") );
    const int nb_edges = edges.shape(0);
    for( int jedge=0; jedge<nb_edges; ++jedge )
    {
      int ip1 = edge_nodes(jedge,0);
      int ip2 = edge_nodes(jedge,1);
      LLPoint2 p1(coords(ip1,LON),coords(ip1,LAT));
      LLPoint2 p2(coords(ip2,LON),coords(ip2,LAT));
      dist(jedge) = arc_in_rad(p1,p2)/DEG_TO_RAD;

      if( p1.lat() == p2.lat() )
      {
        dist_lon_deg(ip1) = dist(jedge);
        dist_lon_deg(ip2) = dist(jedge);
      }
    }
  }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_statistics ( Mesh* mesh) {
  build_statistics(*mesh);
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

