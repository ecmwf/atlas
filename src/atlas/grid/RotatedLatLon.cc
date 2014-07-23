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

ConcreteBuilderT1<Grid,RotatedLatLon> RotatedLatLon_builder;

//void RotatedLatLon::constructFrom(const GridSpec& grid_spec)
//{
//   if (grid_spec.has("Nj")) nptsNS_ = grid_spec.get("Nj");
//   if (grid_spec.has("Ni")) nptsWE_ = grid_spec.get("Ni");
//   if (grid_spec.has("rotated_latitude")) rotated_latitude_ = grid_spec.get("rotated_latitude");
//   if (grid_spec.has("rotated_longitude")) rotated_longitude_ = grid_spec.get("rotated_longitude");
//   if (grid_spec.has("rotated_angle")) rotated_angle_ = grid_spec.get("rotated_angle");
//   if (grid_spec.has("nsIncrement")) nsIncrement_ = grid_spec.get("nsIncrement");
//   if (grid_spec.has("weIncrement")) weIncrement_ = grid_spec.get("weIncrement");
//   if (grid_spec.has("hash"))    hash_ = (std::string)grid_spec.get("hash");
//   grid_spec.get_bounding_box(bbox_);
//   grid_spec.get_points(points_);
//}

RotatedLatLon::RotatedLatLon( const eckit::Params& p )
{
	NOTIMP;
}

RotatedLatLon::~RotatedLatLon()
{
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


void RotatedLatLon::coordinates(std::vector<double>& r ) const
{
	NOTIMP;
}

void RotatedLatLon::coordinates(std::vector<Grid::Point>&) const
{
	NOTIMP;
}

GridSpec* RotatedLatLon::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "RL" << nptsNS_;
   grid_spec->set_short_name(ss.str());
   grid_spec->set("Ni",eckit::Value(nptsWE_));
   grid_spec->set("Nj",eckit::Value(nptsNS_));

   grid_spec->set("rotated_latitude",eckit::Value(rotated_latitude_));
   grid_spec->set("rotated_longitude",eckit::Value(rotated_longitude_));
   grid_spec->set("rotated_angle",eckit::Value(rotated_angle_));
   grid_spec->set("nsIncrement",eckit::Value(nsIncrement_));
   grid_spec->set("weIncrement",eckit::Value(weIncrement_));

   grid_spec->set("hash",eckit::Value(hash_));

   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_points(points_);

   return grid_spec;
}

bool RotatedLatLon::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const RotatedLatLon&>(grid).nptsNS_ != nptsNS_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).nptsWE_ != nptsWE_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).rotated_latitude_ != rotated_latitude_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).rotated_longitude_ != rotated_longitude_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).rotated_angle_ != rotated_angle_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).weIncrement_ != weIncrement_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).hash_ != hash_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const RotatedLatLon&>(grid).points_ != points_) return false;

   return true;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
