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
#include "eckit/value/Value.h"
#include "eckit/geometry/RotateGrid.h"

#include "atlas/GridSpec.h"
#include "atlas/grids/RotatedLatLon.h"

using namespace eckit;
using namespace eckit::geometry;
using namespace std;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,RotatedLatLon,RotatedLatLon::grid_type_str());

RotatedLatLon::RotatedLatLon( const eckit::Params& p )
: south_pole_lat_(0),
  south_pole_lon_(0),
  south_pole_rot_angle_(0),
  nsIncrement_(0),
  weIncrement_(0),
  nptsNS_(0),
  nptsWE_(0)
{
   bbox_ = make_bounding_box(p);

   if( p.has("Nj") )
   {
      nptsNS_ = p["Nj"];
   }

   if( p.has("Ni") )
   {
      nptsWE_ = p["Ni"];
   }

   if( p.has("lat_inc") )
   {
      nsIncrement_ = p["lat_inc"];
   }

   if( p.has("lon_inc") )
   {
      weIncrement_ = p["lon_inc"];
   }

   if( p.has("SouthPoleLat") )
   {
      south_pole_lat_ = p["SouthPoleLat"];
   }

   if( p.has("SouthPoleLon") )
   {
      south_pole_lon_ = p["SouthPoleLon"];
   }

   if( p.has("SouthPoleRotAngle") )
   {
      south_pole_rot_angle_ = p["SouthPoleRotAngle"];
   }

   ASSERT( nptsNS_ > 0 );
   ASSERT( nptsWE_ > 0 );
   ASSERT( nsIncrement_ > 0 );
   ASSERT( weIncrement_ > 0 );
}

RotatedLatLon::~RotatedLatLon()
{
}

std::string RotatedLatLon::shortName() const {

  if( shortName_.empty() )
  {
    std::ostringstream s;
    s <<  grid_type_str()
      << "." << nptsNS_ << "x" << nptsWE_
      << ".R" << south_pole_rot_angle_ << eckit::StrStream::ends;
    shortName_ = s.str();
  }

  return shortName_;
}

void RotatedLatLon::hash(eckit::MD5& md5) const {

  md5.add(grid_type_str());

  md5.add(nptsWE_);
  md5.add(nptsNS_);
  md5.add(south_pole_lat_);
  md5.add(south_pole_lon_);
  md5.add(south_pole_rot_angle_);
  md5.add(nsIncrement_);
  md5.add(weIncrement_);

  bbox_.hash(md5);
}

Grid::Point RotatedLatLon::lonlat(size_t jlon, size_t jlat) const
{
   RotateGrid rotgrid(Grid::Point(south_pole_lon_,south_pole_lat_),south_pole_rot_angle_);

   double plon = bbox_.min().lon(); // west
   double plat = bbox_.max().lat();   // north;
   for( size_t ilat = 0; ilat < nptsNS_; ++ilat) {
      for( size_t ilon = 0; ilon < nptsWE_; ++ilon) {
         if (jlon == ilon && jlat == ilat) {
            return rotgrid.rotate( Grid::Point( plon, plat ) );
         }
         plon += weIncrement_;
      }
      plat += nsIncrement_;
   }
   // should not be here
   return Grid::Point();
}

size_t RotatedLatLon::npts() const
{
   return nptsNS_ * nptsWE_;
}

void RotatedLatLon::lonlat( double pts[] ) const
{
   std::vector<Grid::Point> gpts;
   gpts.resize( npts() );
   lonlat(gpts);

   int c(0);
   for(size_t i = 0; i < gpts.size(); i++) {
      pts[c++] = gpts[i].lon();
      pts[c++] = gpts[i].lat();
   }
}

void RotatedLatLon::lonlat(std::vector<Grid::Point>& points) const
{
   points.resize(npts());
   RotateGrid rotgrid(Grid::Point(south_pole_lon_,south_pole_lat_),south_pole_rot_angle_);

   size_t index = 0;
   double plon = bbox_.min().lon(); // west
   double plat = bbox_.max().lat();   // north;
   for( size_t j = 0; j < nptsNS_; ++j) {
      for( size_t i = 0; i < nptsWE_; ++i) {
         ASSERT( index < points.size() );
         points[index] = rotgrid.rotate( Grid::Point( plon, plat ) );
         index++;
         plon += weIncrement_;
      }
      plat += nsIncrement_;
   }
}

string RotatedLatLon::grid_type() const { return RotatedLatLon::grid_type_str(); }

GridSpec RotatedLatLon::spec() const
{
   GridSpec grid_spec(grid_type());

   grid_spec.set("Ni",eckit::Value(nptsWE_));
   grid_spec.set("Nj",eckit::Value(nptsNS_));

   grid_spec.set("SouthPoleLat",eckit::Value(south_pole_lat_));
   grid_spec.set("SouthPoleLon",eckit::Value(south_pole_lon_));
   grid_spec.set("SouthPoleRotAngle",eckit::Value(south_pole_rot_angle_));
   grid_spec.set("lat_inc",eckit::Value(nsIncrement_));
   grid_spec.set("lon_inc",eckit::Value(weIncrement_));

   grid_spec.set_bounding_box(bbox_);

   return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
