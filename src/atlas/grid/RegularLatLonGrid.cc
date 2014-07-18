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
    /// @todo this function is VERY inneficient -- please rewrite it!

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

double RegularLatLonGrid::computeIncLat() const
{
    /// @note why not simply this??!!
    //       return (bbox_.north() - bbox_.south()) / rows();

    double north = bbox_.top_right().lat();
    double south = bbox_.bottom_left().lat();

    double north_diff_south = 0.0;
    if (north > 0.0 && south > 0.0 ) north_diff_south = north - south;
    else if ( north < 0.0 && south < 0.0) north_diff_south = fabs(north) - fabs(south);
    else north_diff_south  = fabs(north) + fabs(south);

    if (rows() > north_diff_south)
        return north_diff_south/rows();

    // Avoid truncation errors
    long inc_lat = north_diff_south/(rows() + 1) + 0.5;
    return inc_lat;
}

double RegularLatLonGrid::computeIncLon() const
{
    double east  = bbox_.top_right().lon();
    double west  = bbox_.bottom_left().lon();

    if (cols() > (east - west))
        return ((east - west)/cols());

    // Avoid truncation errors
    long inc_lon = ((east - west)/cols() + 0.5 );
    return inc_lon;
}

long RegularLatLonGrid::computeRows(double north, double south, double west, double east) const
{
    if (north > 0.0 && south > 0.0 )
        return (north - south)/incLat() + 1;
    else
        if ( north < 0.0 && south < 0.0)
            return (fabs(north) - fabs(south))/incLat() + 1;

    return (fabs(north) + fabs(south))/incLat() + 1;
}

long RegularLatLonGrid::computeCols(double west, double east) const
{
    return fabs((east - west)/ incLon()) + 1;
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
    if (grid_spec.has("Nj"))      nptsNS_ = grid_spec.get("Nj");
    if (grid_spec.has("Ni"))      nptsWE_ = grid_spec.get("Ni");
    if (grid_spec.has("hash"))    hash_ = (std::string)grid_spec.get("hash");
    grid_spec.get_bounding_box(bbox_);
    grid_spec.get_points(points_);
}

void RegularLatLonGrid::constructFrom(const eckit::Params& p)
{
    nsIncrement_ = p.get("grid_ns");
    weIncrement_ = p.get("grid_ew");

    DEBUG_VAR(nsIncrement_);
    DEBUG_VAR(weIncrement_);

    if( ! p.get("area_s").isNil() )
    {
        bbox_ = BoundBox( p.get("area_n"), p.get("area_s"), p.get("area_e"), p.get("area_w") );
    }

    DEBUG_VAR(bbox_);

    NOTIMP;

}

bool RegularLatLonGrid::compare(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const RegularLatLonGrid&>(grid).nptsNS_ != nptsNS_) return false;
   if ( static_cast<const RegularLatLonGrid&>(grid).nptsWE_ != nptsWE_) return false;
   if ( static_cast<const RegularLatLonGrid&>(grid).hash_ != hash_) return false;
   if ( static_cast<const RegularLatLonGrid&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const RegularLatLonGrid&>(grid).points_ != points_) return false;

   return true;
}

REGISTERIMPL(RegularLatLonGrid,"regular_ll");

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
