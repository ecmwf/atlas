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

#include "eckit/log/Log.h"
#include "eckit/value/Value.h"
#include "atlas/grid/ReducedGaussianGrid.h"
#include "atlas/grid/GridSpec.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ? No.
// NPoints: is this just the grids points, or includes area points, of area does not fit grid
//          assumes it is grid points inclusive of the area.

ReducedGaussianGrid::ReducedGaussianGrid()
: gaussianNumber_(0),nj_(0)
{
}

ReducedGaussianGrid::~ReducedGaussianGrid()
{
//    Log::info() << "Destroy a ReducedGaussianGrid" << std::endl;
}

void ReducedGaussianGrid::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

GridSpec* ReducedGaussianGrid::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   std::stringstream ss; ss << "QG" << gaussianNumber_;
   grid_spec->set_short_name(ss.str());
   grid_spec->set("gaussianNumber",eckit::Value(gaussianNumber_));

   grid_spec->set("hash",eckit::Value(hash_));

   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_rgspec(rgSpec_);
   grid_spec->set_latitudes(latitudes_);
   grid_spec->set_points(points_);

   return grid_spec;
}

void ReducedGaussianGrid::constructFrom(const GridSpec& grid_spec)
{
   if (grid_spec.has("gaussianNumber")) gaussianNumber_ = grid_spec.get("gaussianNumber");
   if (grid_spec.has("hash"))           hash_ = (std::string)grid_spec.get("hash");
   grid_spec.get_bounding_box(bbox_);
   grid_spec.get_rgspec(rgSpec_);
   grid_spec.get_latitudes(latitudes_);
   grid_spec.get_points(points_);
}

bool ReducedGaussianGrid::compare(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const ReducedGaussianGrid&>(grid).gaussianNumber_ != gaussianNumber_) return false;
   if ( static_cast<const ReducedGaussianGrid&>(grid).hash_ != hash_) return false;
   if ( static_cast<const ReducedGaussianGrid&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const ReducedGaussianGrid&>(grid).rgSpec_ != rgSpec_) return false;
   if ( static_cast<const ReducedGaussianGrid&>(grid).latitudes_ != latitudes_) return false;
   if ( static_cast<const ReducedGaussianGrid&>(grid).points_ != points_) return false;

   return true;
}


REGISTERIMPL(ReducedGaussianGrid,"reduced_gg");

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
