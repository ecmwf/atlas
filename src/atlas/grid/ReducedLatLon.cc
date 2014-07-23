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
#include "atlas/grid/ReducedLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,ReducedLatLon> ReducedLatLon_builder;

//void ReducedLatLon::constructFrom(const GridSpec& grid_spec)
//{
//   if (grid_spec.has("Nj")) nptsNS_ = grid_spec.get("Nj");
//   if (grid_spec.has("nsIncrement")) nsIncrement_ =  grid_spec.get("nsIncrement");
//   if (grid_spec.has("hash"))        hash_ = (std::string)grid_spec.get("hash");
//   grid_spec.get_bounding_box(bbox_);
//   grid_spec.get_rgspec(rgSpec_);
//   grid_spec.get_points(points_);
//}

ReducedLatLon::ReducedLatLon( const eckit::Params& p )
{
	NOTIMP;
}

ReducedLatLon::~ReducedLatLon()
{
}

void ReducedLatLon::coordinates(std::vector<double>& r ) const
{
	NOTIMP;
}

void ReducedLatLon::coordinates(std::vector<Grid::Point>&) const
{
	NOTIMP;
}

string ReducedLatLon::gridType() const
{
	return ReducedLatLon::gridTypeStr();
}

GridSpec* ReducedLatLon::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss;

   ss << ReducedLatLon::gridTypeStr() << "_" << nptsNS_;

   grid_spec->set_short_name(ss.str());

   grid_spec->set("Nj",eckit::Value(nptsNS_));
   grid_spec->set("nsIncrement",eckit::Value(nsIncrement_));

   grid_spec->set("hash",eckit::Value(hash_));

   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_rgspec(rgSpec_);
   grid_spec->set_points(points_);

   return grid_spec;
}

bool ReducedLatLon::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const ReducedLatLon&>(grid).nptsNS_ != nptsNS_) return false;
   if ( static_cast<const ReducedLatLon&>(grid).nsIncrement_ != nsIncrement_) return false;
   if ( static_cast<const ReducedLatLon&>(grid).hash_ != hash_) return false;
   if ( static_cast<const ReducedLatLon&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const ReducedLatLon&>(grid).rgSpec_ != rgSpec_) return false;
   if ( static_cast<const ReducedLatLon&>(grid).points_ != points_) return false;

   return true;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
