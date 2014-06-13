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

//#include "eckit/log/Log.h"
#include "atlas/grid/RegularGaussianGrid.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

RegularGaussianGrid::RegularGaussianGrid()
: gaussianNumber_(0),
  the_grid_spec_("regular_gg")
{
//   Log::info() << "RegularGaussianGrid" << std::endl;
}

RegularGaussianGrid::~RegularGaussianGrid()
{
//    Log::info() << "Destroy a RegularGaussianGrid" << std::endl;
}

Grid::Point RegularGaussianGrid::latLon(size_t the_i, size_t the_j) const
{
   long nptsWE = 4 * gaussianNumber_ ;
   long weIncrement = 360.0 / nptsWE;
   for(size_t i = 0 ; i < latitudes_.size(); i++) {

      double plon = bbox_.bottom_left_.lon(); // west_;
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

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
