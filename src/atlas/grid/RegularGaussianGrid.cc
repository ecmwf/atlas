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
#include "atlas/grid/RegularGaussianGrid.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

RegularGaussianGrid::RegularGaussianGrid()
: gaussianNumber_(0),nj_(0)
{
//   Log::info() << "RegularGaussianGrid" << std::endl;
}

RegularGaussianGrid::~RegularGaussianGrid()
{
//    Log::info() << "Destroy a RegularGaussianGrid" << std::endl;
}

Grid::Point RegularGaussianGrid::latLon(size_t the_i, size_t the_j) const
{
    /// @todo this function is VERY inneficient -- please rewrite it!

    long nptsWE = 4 * gaussianNumber_ ;
    long weIncrement = 360.0 / nptsWE;
    for(size_t i = 0 ; i < latitudes_.size(); i++) {

        double plon = bbox_.bottom_left().lon(); // west_;
        for( size_t j = 0; j < nptsWE; ++j) {
            if ( i== the_i && j== the_j) {
                return  Point( latitudes_[i], plon );
            }
            plon += weIncrement;
        }
    }

    return Grid::Point();
}

void RegularGaussianGrid::coordinates( Grid::Coords& r ) const
{
   ASSERT( r.size() == points_.size() );

   for( size_t i = 0; i < points_.size(); ++i )
   {
      r.lat(i) = points_[i].lat();
      r.lon(i) = points_[i].lon();
   }
}

GridSpec* RegularGaussianGrid::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "GG" << gaussianNumber_;
   grid_spec->set_short_name(ss.str());
   grid_spec->set("Nj",eckit::Value(nj_));
   grid_spec->set("gaussianNumber",eckit::Value(gaussianNumber_));

   grid_spec->set("hash",eckit::Value(hash_));
   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_latitudes(latitudes_);
   grid_spec->set_points(points_);

   return grid_spec;
}

void RegularGaussianGrid::constructFrom(const GridSpec& grid_spec)
{
   if (grid_spec.has("Nj"))           nj_ = grid_spec.get("Nj");
   if (grid_spec.has("gaussianNumber")) gaussianNumber_ = grid_spec.get("gaussianNumber");
   if (grid_spec.has("hash"))           hash_ = (std::string)grid_spec.get("hash");
   grid_spec.get_bounding_box(bbox_);
   grid_spec.get_latitudes(latitudes_);
   grid_spec.get_points(points_);
}

bool RegularGaussianGrid::compare(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const RegularGaussianGrid&>(grid).gaussianNumber_ != gaussianNumber_) return false;
   if ( static_cast<const RegularGaussianGrid&>(grid).nj_ != nj_) return false;
   if ( static_cast<const RegularGaussianGrid&>(grid).hash_ != hash_) return false;
   if ( static_cast<const RegularGaussianGrid&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const RegularGaussianGrid&>(grid).latitudes_ != latitudes_) return false;
   if ( static_cast<const RegularGaussianGrid&>(grid).points_ != points_) return false;

   return true;
}

REGISTERIMPL(RegularGaussianGrid,"regular_gg");

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
