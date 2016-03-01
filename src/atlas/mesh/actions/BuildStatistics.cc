/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <set>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/filesystem/PathName.h"
#include "atlas/atlas_config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/Constants.h"
#include "atlas/array/ArrayView.h"
#include "atlas/util/runtime/ErrorHandling.h"
#include "atlas/util/parallel/mpl/Checksum.h"

using eckit::geometry::LLPoint2;
using eckit::geometry::lonlat_to_3d;
namespace atlas {
namespace mesh {
namespace actions {

static const double DEG_TO_RAD = M_PI/180.;

namespace {

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

double quad_quality( const LLPoint2& p1, const LLPoint2& p2, const LLPoint2& p3, const LLPoint2& p4)
{
  // see http://geuz.org/gmsh/doc/preprints/gmsh_quad_preprint.pdf

  double xyz1[3];
  double xyz2[3];
  double xyz3[3];
  double xyz4[3];

  lonlat_to_3d(p1.lon(),p1.lat(),xyz1,1.,0.);
  lonlat_to_3d(p2.lon(),p2.lat(),xyz2,1.,0.);
  lonlat_to_3d(p3.lon(),p3.lat(),xyz3,1.,0.);
  lonlat_to_3d(p4.lon(),p4.lat(),xyz4,1.,0.);

  double l2m1[3];
  l2m1[0] = xyz2[0]-xyz1[0];
  l2m1[2] = xyz2[1]-xyz1[1];
  l2m1[1] = xyz2[2]-xyz1[2];

  double l3m2[3];
  l3m2[0] = xyz3[0]-xyz2[0];
  l3m2[2] = xyz3[1]-xyz2[1];
  l3m2[1] = xyz3[2]-xyz2[2];

  double l4m3[3];
  l4m3[0] = xyz4[0]-xyz3[0];
  l4m3[2] = xyz4[1]-xyz3[1];
  l4m3[1] = xyz4[2]-xyz3[2];

  double l1m4[3];
  l1m4[0] = xyz1[0]-xyz4[0];
  l1m4[2] = xyz1[1]-xyz4[1];
  l1m4[1] = xyz1[2]-xyz4[2];

  double norm_l2m1 = std::sqrt( l2m1[0]*l2m1[0] + l2m1[1]*l2m1[1] + l2m1[2]*l2m1[2] );
  double norm_l3m2 = std::sqrt( l3m2[0]*l3m2[0] + l3m2[1]*l3m2[1] + l3m2[2]*l3m2[2] );
  double norm_l4m3 = std::sqrt( l4m3[0]*l4m3[0] + l4m3[1]*l4m3[1] + l4m3[2]*l4m3[2] );
  double norm_l1m4 = std::sqrt( l1m4[0]*l1m4[0] + l1m4[1]*l1m4[1] + l1m4[2]*l1m4[2] );

  double dot_l4m1_l2m1 = -l1m4[0]*l2m1[0]-l1m4[1]*l2m1[1]-l1m4[2]*l2m1[2];
  double dot_l1m2_l3m2 = -l2m1[0]*l3m2[0]-l2m1[1]*l3m2[1]-l2m1[2]*l3m2[2];
  double dot_l2m3_l4m3 = -l3m2[0]*l4m3[0]-l3m2[1]*l4m3[1]-l3m2[2]*l4m3[2];
  double dot_l3m4_l1m4 = -l4m3[0]*l1m4[0]-l4m3[1]*l1m4[1]-l4m3[2]*l1m4[2];

  // Angles at each quad corner
  double a1 = std::acos( dot_l4m1_l2m1 / ( norm_l1m4 * norm_l2m1 ) );
  double a2 = std::acos( dot_l1m2_l3m2 / ( norm_l2m1 * norm_l3m2 ) );
  double a3 = std::acos( dot_l2m3_l4m3 / ( norm_l3m2 * norm_l4m3 ) );
  double a4 = std::acos( dot_l3m4_l1m4 / ( norm_l4m3 * norm_l1m4 ) );

  double max_inner = std::max( std::max( std::max(
    std::abs(M_PI_2-a1), std::abs(M_PI_2-a2) ), std::abs(M_PI_2-a3) ), std::abs(M_PI_2-a4));

  return std::max(1. - M_2_PI*max_inner, 0.);
}

}

void build_statistics( Mesh& mesh )
{
  const double radius_km = util::Earth::radiusInMeters()*1e-3;

  mesh::Nodes& nodes = mesh.nodes();
  array::ArrayView<double,2> lonlat ( nodes.lonlat() );

  if( mesh.edges().size() )
  {
    if( ! mesh.edges().has_field("arc_length") )
      mesh.edges().add( field::Field::create<double>("arc_length",array::make_shape(mesh.edges().size())) );
    array::ArrayView<double,1> dist ( mesh.edges().field("arc_length") );
    const mesh::HybridElements::Connectivity &edge_nodes = mesh.edges().node_connectivity();

    const int nb_edges = mesh.edges().size();
    for( int jedge=0; jedge<nb_edges; ++jedge )
    {
      int ip1 = edge_nodes(jedge,0);
      int ip2 = edge_nodes(jedge,1);
      LLPoint2 p1(lonlat(ip1,internals::LON),lonlat(ip1,internals::LAT));
      LLPoint2 p2(lonlat(ip2,internals::LON),lonlat(ip2,internals::LAT));
      dist(jedge) = arc_in_rad(p1,p2)*radius_km;
    }
  }

  std::ofstream ofs;
  eckit::PathName stats_path("stats.txt");
  int idt = 10;
  if( eckit::mpi::size() == 1 )
  {
    ofs.open( stats_path.localPath(), std::ofstream::out );
    ofs << "# STATISTICS rho (min_length/max_length), eta (quality) \n";
    ofs << std::setw(idt) << "# rho";
    ofs << std::setw(idt) << "eta";
    ofs << "\n";
    ofs.close();
  }

  // Cell statistics
  {
    if( eckit::mpi::size() == 1 )
      ofs.open( stats_path.localPath(), std::ofstream::app );

    array::ArrayView<double,1> rho ( mesh.cells().add( field::Field::create<double>("stats_rho",array::make_shape(mesh.cells().size()) ) ) );
    array::ArrayView<double,1> eta ( mesh.cells().add( field::Field::create<double>("stats_eta",array::make_shape(mesh.cells().size()) ) ) );

    for( size_t jtype=0; jtype<mesh.cells().nb_types(); ++jtype )
    {
      const mesh::Elements& elements = mesh.cells().elements(jtype);
      const mesh::Elements::Connectivity& elem_nodes = elements.node_connectivity();
      const size_t nb_elems = elements.size();

      if( elements.element_type().name() == "Triangle" )
      {
        for( size_t jelem=0; jelem<nb_elems; ++jelem )
        {
          size_t ielem = elements.begin() + jelem;
          size_t ip1 = elem_nodes(jelem,0);
          size_t ip2 = elem_nodes(jelem,1);
          size_t ip3 = elem_nodes(jelem,2);
          LLPoint2 p1(lonlat(ip1,internals::LON),lonlat(ip1,internals::LAT));
          LLPoint2 p2(lonlat(ip2,internals::LON),lonlat(ip2,internals::LAT));
          LLPoint2 p3(lonlat(ip3,internals::LON),lonlat(ip3,internals::LAT));

          double l12 = arc_in_rad(p1,p2)/DEG_TO_RAD;
          double l23 = arc_in_rad(p2,p3)/DEG_TO_RAD;
          double l31 = arc_in_rad(p3,p1)/DEG_TO_RAD;

          double min_length = std::min(std::min(l12,l23),l31);
          double max_length = std::max(std::max(l12,l23),l31);
          rho(ielem) = min_length/max_length;

          double s = 0.5*(l12 + l23 + l31);
          double area = std::sqrt(s * (s - l12) * (s - l23) * (s - l31));

          // see http://www.gidhome.com/component/manual/referencemanual/preprocessing/mesh_menu/mesh_quality
          eta(ielem) = (4*area*std::sqrt(3.))/( std::pow(l12,2)+std::pow(l23,2)+std::pow(l31,2) );

          if( eckit::mpi::size() == 1 )
          {
            ofs << std::setw(idt) << rho[ielem]
                << std::setw(idt) << eta[ielem]
                << "\n";
          }
        }
      }
      if( elements.element_type().name() == "Quadrilateral" )
      {
        for( size_t jelem=0; jelem<nb_elems; ++jelem )
        {
          size_t ielem = elements.begin() + jelem;
          size_t ip1 = elem_nodes(jelem,0);
          size_t ip2 = elem_nodes(jelem,1);
          size_t ip3 = elem_nodes(jelem,2);
          size_t ip4 = elem_nodes(jelem,3);

          LLPoint2 p1(lonlat(ip1,internals::LON),lonlat(ip1,internals::LAT));
          LLPoint2 p2(lonlat(ip2,internals::LON),lonlat(ip2,internals::LAT));
          LLPoint2 p3(lonlat(ip3,internals::LON),lonlat(ip3,internals::LAT));
          LLPoint2 p4(lonlat(ip4,internals::LON),lonlat(ip4,internals::LAT));

          eta(ielem) = quad_quality(p1,p2,p3,p4);

          double l12 = arc_in_rad(p1,p2);
          double l23 = arc_in_rad(p2,p3);
          double l34 = arc_in_rad(p3,p4);
          double l41 = arc_in_rad(p4,p1);

          double min_length = std::min(std::min(std::min(l12,l23),l34),l41);
          double max_length = std::max(std::max(std::max(l12,l23),l34),l41);
          rho(ielem) = min_length/max_length;

          if( eckit::mpi::size() == 1 )
          {
            ofs << std::setw(idt) << rho[ielem]
                << std::setw(idt) << eta[ielem]
                << "\n";
          }
        }

      }
    }
    if( eckit::mpi::size() == 1 )
      ofs.close();
  }

  eckit::PathName dual_stats_path("dual_stats.txt");
  if( eckit::mpi::size() == 1 )
  {
    ofs.open( dual_stats_path.localPath(), std::ofstream::out );
    ofs << "# STATISTICS dual_area \n";
    ofs << std::setw(idt) << "# area";
    ofs << "\n";
  }

  if( nodes.has_field("dual_volumes") )
  {
    array::ArrayView<double,1> dual_volumes ( nodes.field("dual_volumes") );
    array::ArrayView<double,1> dual_delta_sph  ( nodes.add( field::Field::create<double>( "dual_delta_sph", array::make_shape(nodes.size(),1) ) ) );

    for( size_t jnode=0; jnode<nodes.size(); ++jnode )
    {
      const double lat = lonlat(jnode,internals::LAT)*DEG_TO_RAD;
      const double hx = radius_km*std::cos(lat)*DEG_TO_RAD;
      const double hy = radius_km*DEG_TO_RAD;
      dual_delta_sph(jnode) = std::sqrt(dual_volumes(jnode)*hx*hy);
    }

    if( eckit::mpi::size() == 1 )
    {
      for( size_t jnode=0; jnode<nodes.size(); ++jnode )
      {
        ofs << std::setw(idt) << dual_delta_sph(jnode)
            << "\n";
      }
    }
  }
  if( eckit::mpi::size() == 1 )
    ofs.close();


}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_statistics ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_statistics(*mesh) );
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

