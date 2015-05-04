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

#include "atlas/GridSpec.h"
#include "atlas/grids/PolarStereoGraphic.h"

using namespace eckit;
using namespace eckit::geometry;
using namespace std;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,PolarStereoGraphic,PolarStereoGraphic::grid_type_str());

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
   first_grid_pt_.assign(lon,lat);


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
   GridSpec grid_spec(grid_type());

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
   grid_spec.set("resolutionAndComponentFlag",eckit::Value(resolutionAndComponentFlag_));

   // Bounding box can be computed
   // grid_spec.set_bounding_box(bounding_box());

   return grid_spec;
}

std::string PolarStereoGraphic::shortName() const {

  if( shortName_.empty() )
  {
    std::ostringstream s;
    s <<  grid_type_str()
      << "." << npts_xaxis_ << "x" << npts_yaxis_
      << eckit::StrStream::ends;
    shortName_ = s.str();
  }

  return shortName_;
}

Grid::uid_t PolarStereoGraphic::unique_id() const {

  if (uid_.empty()) {
    eckit::StrStream os;
    os << shortName() << "." << hash() << eckit::StrStream::ends;
    uid_ = std::string(os);
  }

  return uid_;
}

MD5::digest_t PolarStereoGraphic::hash() const {

  if (hash_.empty()) {

    eckit::MD5 md5;

    md5.add(npts_xaxis_);
    md5.add(npts_yaxis_);

    md5.add(x_grid_length_);
    md5.add(y_grid_length_);

    md5.add(first_grid_pt_.lat());
    md5.add(first_grid_pt_.lon());


    md5.add(orientationOfTheGrid_);
    md5.add(southPoleOnProjectionPlane_);

    md5.add(earth_is_oblate_);
    md5.add(semi_major_);
    md5.add(semi_minor_);
    md5.add(radius_);

    // "LaD" ???
    // "resolutionAndComponentFlag" ???

    hash_ = md5.digest();
  }

  return hash_;
}

size_t PolarStereoGraphic::npts() const
{
   return npts_xaxis_ * npts_yaxis_;
}

void PolarStereoGraphic::lonlat( double pts[] ) const
{
   std::vector<Grid::Point> gpts;
   lonlat(gpts);

   int c(0);
   for(size_t i = 0; i < gpts.size(); i++) {
      pts[c++] = gpts[i].lon();
      pts[c++] = gpts[i].lat();
   }
}

//Enable to test against grib iterator
//#define GRIB_COMPAT 1

void PolarStereoGraphic::lonlat(std::vector<Grid::Point>& points) const
{
   points.resize(npts());

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
   BoundBox bbox = bounding_box();
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

BoundBox PolarStereoGraphic::bounding_box() const
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

string PolarStereoGraphic::grid_type() const
{
   return PolarStereoGraphic::grid_type_str();
}

//-----------------------------------------------------------------------------


} // namespace grids
} // namespace atlas
