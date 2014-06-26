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

//#include "eckit/log/Log.h"
#include "eckit/value/Value.h"

#include "atlas/grid/RegularLatLonGrid.h"
#include "atlas/grid/GridSpec.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

RegularLatLonGrid::RegularLatLonGrid()
:  nsIncrement_(0),
   weIncrement_(0),
   nptsNS_(0),
   nptsWE_(0)
{
//   Log::info() << "RegularLatLonGrid" << std::endl;
}

RegularLatLonGrid::~RegularLatLonGrid()
{
//    Log::info() << "Destroy a RegularLatLonGrid" << std::endl;
}

Grid::Point RegularLatLonGrid::latLon(size_t the_i, size_t the_j) const
{
   double plon = bbox_.bottom_left_.lon(); // west
   double plat = bbox_.top_right_.lat();   // north;
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


void RegularLatLonGrid::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

GridSpec* RegularLatLonGrid::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "LL" << nptsNS_ << "_" << nptsWE_;
   grid_spec->set_short_name(ss.str());

   grid_spec->set("Nj",eckit::Value(nptsNS_));
   grid_spec->set("Ni",eckit::Value(nptsWE_));

   grid_spec->set("hash",eckit::Value(hash_));
   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_points(points_);

   return grid_spec;
}

void RegularLatLonGrid::constructFrom(const GridSpec& grid_spec)
{
   nptsNS_ = grid_spec.get("Nj");
   nptsWE_ = grid_spec.get("Ni");
   hash_ = (std::string) grid_spec.get("hash");
   grid_spec.get_bounding_box(bbox_);
   grid_spec.get_points(points_);
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
