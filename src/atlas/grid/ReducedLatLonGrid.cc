/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

// ==================================================================================
// gribs use the following convention: (from Shahram)
//
// Horizontally:  Points scan in the +i (+x) direction
// Vertically:    Points scan in the -j (-y) direction
//
// The way I verified this was to look at our SAMPLE files (which IFS uses).
// I also verified that IFS does not modify the scanning modes
// so whatever the samples say, is the convention
// ==================================================================================

#include <stdexcept>

#include "eckit/value/Value.h"
//#include "eckit/log/Log.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/ReducedLatLonGrid.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

ReducedLatLonGrid::ReducedLatLonGrid()
:  nsIncrement_(0),
   nptsNS_(0)
{
//   Log::info() << "ReducedLatLonGrid" << std::endl;
}

ReducedLatLonGrid::~ReducedLatLonGrid()
{
//    Log::info() << "Destroy a ReducedLatLonGrid" << std::endl;
}

void ReducedLatLonGrid::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

GridSpec* ReducedLatLonGrid::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "RedLL" << nptsNS_;
   grid_spec->set_short_name(ss.str());
   grid_spec->set("Nj",eckit::Value(nptsNS_));
   grid_spec->set("nsIncrement",eckit::Value(nsIncrement_));

   grid_spec->set("hash",eckit::Value(hash_));

   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_rgspec(rgSpec_);
   grid_spec->set_points(points_);

   return grid_spec;
}

void ReducedLatLonGrid::constructFrom(const GridSpec& grid_spec)
{
   if (grid_spec.has("Nj")) nptsNS_ = grid_spec.get("Nj");
   if (grid_spec.has("nsIncrement")) nsIncrement_ =  grid_spec.get("nsIncrement");
   if (grid_spec.has("hash"))        hash_ = (std::string)grid_spec.get("hash");
   grid_spec.get_bounding_box(bbox_);
   grid_spec.get_rgspec(rgSpec_);
   grid_spec.get_points(points_);
}

bool ReducedLatLonGrid::compare(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const ReducedLatLonGrid&>(grid).nptsNS_ != nptsNS_) return false;
   if ( static_cast<const ReducedLatLonGrid&>(grid).nsIncrement_ != nsIncrement_) return false;
   if ( static_cast<const ReducedLatLonGrid&>(grid).hash_ != hash_) return false;
   if ( static_cast<const ReducedLatLonGrid&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const ReducedLatLonGrid&>(grid).rgSpec_ != rgSpec_) return false;
   if ( static_cast<const ReducedLatLonGrid&>(grid).points_ != points_) return false;

   return true;
}

REGISTERIMPL(ReducedLatLonGrid,"reduced_ll");

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
