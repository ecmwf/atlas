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
: iScansPositively_(true),
  jScansPositively_(false),
  npts_xaxis_(0),
  npts_yaxis_(0),
  x_grid_length_(0),
  y_grid_length_(0),
  resolutionAndComponentFlag_(8), // default 8 assumes earth is spherical
  orientationOfTheGrid_(0),
  lad_(60),
  southPoleOnProjectionPlane_(false),
  earth_is_oblate_(false),
  radius_(6371229),
  semi_major_(6378137),
  semi_minor_(6356752.3),
  e_(0.081819191)
{
   if( !p.get("hash").isNil() )
      hash_ = p["hash"].as<std::string>();

   if( p.has("iScansPositively") )
   {
      iScansPositively_ = p["iScansPositively"];
   }

   if( p.has("jScansPositively") )
   {
      jScansPositively_ = p["jScansPositively"];
   }

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

   if( p.has("orientationOfTheGrid") )
   {
      orientationOfTheGrid_ = p["orientationOfTheGrid"];
   }

   if( p.has("LaD") )
   {
      lad_ = p["LaD"];
   }

   // resolutionAndComponentFlag is used to indicate sphere/oblate
   // But this is extracted separately, to avoid having to mess with bits.
   // Stored since needed to match geography hash, when gridSpec written to grib
   if( p.has("resolutionAndComponentFlag") )
   {
      resolutionAndComponentFlag_ = p["resolutionAndComponentFlag"];
   }

   double lat = 0;
   double lon = 0;
   if( p.has("latitudeOfFirstGridPoint") )
   {
      lat =  p["latitudeOfFirstGridPoint"];
   }
   if( p.has("longitudeOfFirstGridPoint") )
   {
      lon =  p["longitudeOfFirstGridPoint"];
   }
   first_grid_pt_.assign(lat,lon);


   if (p.has("southPoleOnProjectionPlane")) {
      southPoleOnProjectionPlane_ = p["southPoleOnProjectionPlane"];
   }

   if (p.has("earthIsOblate")) {
      earth_is_oblate_ = p["earthIsOblate"];
   }


   if (earth_is_oblate_) {

      if (p.has("earthMajorAxis")) {
         semi_major_ = p["earthMajorAxis"];
      }

      if (p.has("earthMinorAxis")) {
         semi_minor_ = p["earthMinorAxis"];
      }

      e_ = sqrt( 1.0 - ((double) semi_minor_*semi_minor_) / ((double) semi_major_*semi_major_));
   }
   else {
      if (p.has("radius")) {
         radius_ = p["radius"];
      }
   }

   ASSERT( npts_xaxis_ > 0);
   ASSERT( npts_yaxis_ > 0);
   ASSERT( x_grid_length_ > 0);
   ASSERT( y_grid_length_ > 0);
   ASSERT( semi_major_ > semi_minor_);
   ASSERT( e_ < 1.0);
   ASSERT( lad_ > 0);
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
   grid_spec.set("latitudeOfFirstGridPoint",eckit::Value(first_grid_pt_.lat()));
   grid_spec.set("longitudeOfFirstGridPoint",eckit::Value(first_grid_pt_.lon()));
   grid_spec.set("iScansPositively",eckit::Value(iScansPositively_));
   grid_spec.set("jScansPositively",eckit::Value(jScansPositively_));

   grid_spec.set("orientationOfTheGrid",eckit::Value(orientationOfTheGrid_));
   grid_spec.set("LaD",eckit::Value(lad_));
   grid_spec.set("southPoleOnProjectionPlane",eckit::Value(southPoleOnProjectionPlane_));
   grid_spec.set("earthIsOblate",eckit::Value(earth_is_oblate_));
   grid_spec.set("radius",eckit::Value(radius_));
   grid_spec.set("earthMajorAxis",eckit::Value(semi_major_));
   grid_spec.set("earthMinorAxis",eckit::Value(semi_minor_));
   grid_spec.set("numberOfDataPoints",eckit::Value(nPoints()));
   grid_spec.set("resolutionAndComponentFlag",eckit::Value(resolutionAndComponentFlag_));

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

//Enable to test against grib iterator
//#define GRIB_COMPAT 1

void PolarStereoGraphic::coordinates(std::vector<Grid::Point>& points) const
{
   ASSERT( points.size() == nPoints() );

   PolarStereoGraphicProj ps(southPoleOnProjectionPlane_,earth_is_oblate_,orientationOfTheGrid_);
   if (earth_is_oblate_) {
      ps.set_radius(semi_major_);
      ps.set_eccentricity(e_);
   }
   else {
      ps.set_radius(radius_);
   }

#ifndef GRIB_COMPAT
   // Points go from North,West --> South,east
   Grid::BoundBox bbox = boundingBox();
   const double north = bbox.north();
   const double west  = bbox.west();

   Point2 first_pt_on_plane = ps.map_to_plane(eckit::geometry::LLPoint2(north,west));
   double x = first_pt_on_plane[0];
   double y = first_pt_on_plane[1];
   size_t k = 0;
   for(size_t j = 0; j < npts_yaxis_; j++) {
      x = first_pt_on_plane[0];
      for(size_t i = 0; i < npts_xaxis_; i++) {

         points[k++] = ps.map_to_spherical(x,y);
         x += x_grid_length_;
      }
      y -= y_grid_length_;
   }
#else

   long x_grid_length = (iScansPositively_) ? x_grid_length_ : -x_grid_length_;
   long y_grid_length = (jScansPositively_) ? y_grid_length_ : -y_grid_length_;

   Point2 first_pt_on_plane = ps.map_to_plane(first_grid_pt_);
   double x = first_pt_on_plane[0];
   double y = first_pt_on_plane[1];
   size_t k = 0;
   for(size_t j = 0; j < npts_yaxis_; j++) {
      x = first_pt_on_plane[0];
      for(size_t i = 0; i < npts_xaxis_; i++) {

         points[k++] = ps.map_to_spherical(x,y);
         x += x_grid_length;
      }
      y += y_grid_length;
   }
#endif
}

Grid::BoundBox PolarStereoGraphic::boundingBox() const
{
   // Map first_grid_pt_ to the plane.
   PolarStereoGraphicProj ps(southPoleOnProjectionPlane_,earth_is_oblate_,orientationOfTheGrid_);
   if (earth_is_oblate_) {
      ps.set_radius(semi_major_);
      ps.set_eccentricity(e_);
   }
   else {
      ps.set_radius(radius_);
   }

   Point2 first_pt_on_plane = ps.map_to_plane(first_grid_pt_);

   // *Depending* on the scanning mode, find the last point on the plane
   // Note: without the scanning mode we can not find the correct bounding box.
   double last_point_x = first_pt_on_plane[0] + npts_xaxis_ * x_grid_length_;
   double last_point_y = first_pt_on_plane[1] + npts_yaxis_ * y_grid_length_;

   last_point_x = (iScansPositively_) ? last_point_x : -last_point_x;
   last_point_y = (jScansPositively_) ? last_point_y : -last_point_y;

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
