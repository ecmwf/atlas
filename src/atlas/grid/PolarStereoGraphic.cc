/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/value/Value.h"
#include "eckit/geometry/PolarStereoGraphicProj.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/PolarStereoGraphic.h"

using namespace eckit;
using namespace eckit::geometry;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,PolarStereoGraphic> PolarStereoGraphic_builder( PolarStereoGraphic::gridTypeStr() );

PolarStereoGraphic::PolarStereoGraphic( const eckit::Params& p )
: npts_xaxis_(0),
  npts_yaxis_(0),
  x_grid_length_(0),
  y_grid_length_(0),
  lov_(0),
  lad_(60),
  north_pole_on_projection_plane_(true),
  spherical_earth_(true)
{
   if( !p.get("hash").isNil() )
      hash_ = p["hash"].as<std::string>();

   if( p.has("Nx") )
   {
      npts_xaxis_ = p["Nx"];
   }

   if( p.has("Ny") )
   {
      npts_yaxis_ = p["Ny"];
   }

   if( p.has("Dx") )
   {
      x_grid_length_ = p["Dx"];
   }

   if( p.has("Dy") )
   {
      y_grid_length_ = p["Dy"];
   }

   if( p.has("Lov") )
   {
      lov_ = p["Lov"];
   }

   if( p.has("Lad") )
   {
      lad_ = p["Lad"];
   }

   double lat = 0;
   double lon = 0;
   if( p.has("La1") )
   {
      lat =  p["La1"];
   }
   if( p.has("Lo1") )
   {
      lon =  p["Lo1"];
   }
   first_grid_pt_.assign(lat,lon);


   if (p.has("north_pole_on_projection_plane")) {
      north_pole_on_projection_plane_ = p["north_pole_on_projection_plane"];
   }

   if (p.has("spherical_earth")) {
      spherical_earth_ = p["spherical_earth"];
   }

   ASSERT( npts_xaxis_ > 0);
   ASSERT( npts_yaxis_ > 0);
   ASSERT( x_grid_length_ > 0);
   ASSERT( y_grid_length_ > 0);
   ASSERT( lov_ > 0);
   ASSERT( lad_ > 0);

   // North Pole projection, cant project point on the south pole, and vice versa
   RealCompare<double> isEqual(degrees_eps());
   ASSERT( north_pole_on_projection_plane_ && !(isEqual(lat,-90.0) && isEqual(lon,0)));
   ASSERT( !north_pole_on_projection_plane_ && !(isEqual(lat,90.0) && isEqual(lon,0)));

   // Bounding box is computed
   // bbox_ = makeBBox(p);
}

PolarStereoGraphic::~PolarStereoGraphic()
{
}

GridSpec PolarStereoGraphic::spec() const
{
   GridSpec grid_spec(gridType());

   grid_spec.uid(uid());
   grid_spec.set("Nx",eckit::Value(npts_xaxis_));
   grid_spec.set("Ny",eckit::Value(npts_yaxis_));

   grid_spec.set("Dx",eckit::Value(x_grid_length_));
   grid_spec.set("Dy",eckit::Value(y_grid_length_));
   grid_spec.set("Lov",eckit::Value(lov_));
   grid_spec.set("La1",eckit::Value(first_grid_pt_.lat()));
   grid_spec.set("Lo1",eckit::Value(first_grid_pt_.lon()));
   grid_spec.set("Lad",eckit::Value(lad_));
   grid_spec.set("north_pole_on_projection_plane",eckit::Value(north_pole_on_projection_plane_));
   grid_spec.set("spherical_earth",eckit::Value(spherical_earth_));

   grid_spec.set("hash",eckit::Value(hash_));

   // Bounding box can be computed
   // grid_spec.set_bounding_box(boundingBox());

   return grid_spec;
}


string PolarStereoGraphic::uid() const
{
   std::stringstream ss;
   ss << gridTypeStr() << "_" << npts_xaxis_ << "_" << npts_yaxis_;
   return ss.str();
}

size_t PolarStereoGraphic::nPoints() const
{
   return npts_xaxis_ * npts_yaxis_;
}

void PolarStereoGraphic::coordinates(std::vector<double>& pts ) const
{
   ASSERT( pts.size() && pts.size()%2 == 0 );

   std::vector<Grid::Point> gpts;
   gpts.resize( pts.size()/ 2);
   coordinates(gpts);

   for(size_t i = 0; i < gpts.size(); i++) {
      pts[i] = gpts[i].lat();
      pts[i+1] = gpts[i].lon();
   }
}

void PolarStereoGraphic::coordinates(std::vector<Grid::Point>& points) const
{
   ASSERT( points.size() == nPoints() );

   PolarStereoGraphicProj ps(north_pole_on_projection_plane_,spherical_earth_,lov_);

   Point2 first_pt_on_plane = ps.map_to_plane(first_grid_pt_);
   double x = first_pt_on_plane[0];
   double y = first_pt_on_plane[1];
   size_t k = 0;
   for(size_t j = 0; j < npts_yaxis_; j++) {
      for(size_t i = 0; i < npts_xaxis_; i++) {

         points[k] = ps.map_to_spherical(x,y);
         k++;
         x += x_grid_length_;
      }
      y += y_grid_length_;
   }
}

Grid::BoundBox PolarStereoGraphic::boundingBox() const
{
   // One point first_grid_pt_
   // Map this to the plane.
   PolarStereoGraphicProj ps(north_pole_on_projection_plane_,spherical_earth_,lov_);
   Point2 first_pt_on_plane = ps.map_to_plane(first_grid_pt_);

   // Find the last point on the plane
   double last_point_x = first_pt_on_plane[0] + npts_xaxis_ * x_grid_length_;
   double last_point_y = first_pt_on_plane[1] + npts_yaxis_ * y_grid_length_;

   // Transform this back into spherical co-ordinates
   Point last_pt_in_sp = ps.map_to_spherical(last_point_x,last_point_y);

   double top = std::max(first_grid_pt_.lat(),last_pt_in_sp.lat());
   double bottom = std::min(first_grid_pt_.lat(),last_pt_in_sp.lat());
   double right = std::max(first_grid_pt_.lon(),last_pt_in_sp.lon());
   double left = std::min(first_grid_pt_.lon(),last_pt_in_sp.lon());

   return BoundBox(top,bottom,right,left);
}

string PolarStereoGraphic::gridType() const
{
   return PolarStereoGraphic::gridTypeStr();
}

bool PolarStereoGraphic::same(const Grid& grid) const
{
   return spec() == grid.spec();
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
