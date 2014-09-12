/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/value/Value.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/RotatedLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,RotatedLatLon> RotatedLatLon_builder( RotatedLatLon::gridTypeStr() );

RotatedLatLon::RotatedLatLon( const eckit::Params& p )
: south_pole_lat_(0),
  south_pole_lon_(0),
  south_pole_rot_angle_(0),
  nsIncrement_(0),
  weIncrement_(0),
  nptsNS_(0),
  nptsWE_(0)
{
   if( !p.get("hash").isNil() )
      hash_ = p["hash"].as<std::string>();

   bbox_ = makeBBox(p);

   if( p.has("Nj") )
   {
      nptsNS_ = p["Nj"];
   }

   if( p.has("Ni") )
   {
      nptsWE_ = p["Ni"];
   }

   if( p.has("grid_lat_inc") )
   {
      nsIncrement_ = p["grid_lat_inc"];
   }

   if( p.has("grid_lon_inc") )
   {
      weIncrement_ = p["grid_lon_inc"];
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

string RotatedLatLon::uid() const
{
	std::stringstream ss;
	ss << gridTypeStr() << "_" << nptsNS_;
	return ss.str();
}

Grid::Point RotatedLatLon::latLon(size_t the_i, size_t the_j) const
{
   double plon = bbox_.bottom_left().lon(); // west
   double plat = bbox_.top_right().lat();   // north;
   for( size_t j = 0; j <= nptsNS_; ++j) {
      for( size_t i = 0; i <= nptsWE_; ++i) {
         if (the_i == i && the_j == j) {
            return Grid::Point( plat, plon );
         }
         plon += weIncrement_;
      }
      plat += nsIncrement_;
   }
   return Grid::Point();
}


size_t RotatedLatLon::nPoints() const
{
   // ?? TODO do we take bbox_ in to account
   return nptsNS_ * nptsWE_;
}


void RotatedLatLon::coordinates(std::vector<double>& pts ) const
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

void RotatedLatLon::coordinates(std::vector<Grid::Point>&) const
{
	NOTIMP;
}

string RotatedLatLon::gridType() const
{
	return RotatedLatLon::gridTypeStr();
}

GridSpec RotatedLatLon::spec() const
{
   GridSpec grid_spec(gridType());

   grid_spec.uid(uid());
   grid_spec.set("Ni",eckit::Value(nptsWE_));
   grid_spec.set("Nj",eckit::Value(nptsNS_));

   grid_spec.set("SouthPoleLat",eckit::Value(south_pole_lat_));
   grid_spec.set("SouthPoleLon",eckit::Value(south_pole_lon_));
   grid_spec.set("SouthPoleRotAngle",eckit::Value(south_pole_rot_angle_));
   grid_spec.set("grid_lat_inc",eckit::Value(nsIncrement_));
   grid_spec.set("grid_lon_inc",eckit::Value(weIncrement_));

   grid_spec.set("hash",eckit::Value(hash_));

   grid_spec.set_bounding_box(bbox_);

   return grid_spec;
}

bool RotatedLatLon::same(const Grid& grid) const
{
	return spec() == grid.spec();
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
