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
#include "atlas/grid/RegularLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,RegularLatLon> RegularLatLon_builder;


//void RegularLatLon::RegularLatLon(const GridSpec& grid_spec)
//{
//    if (grid_spec.has("Nj"))      nptsNS_ = grid_spec.get("Nj");
//    if (grid_spec.has("Ni"))      nptsWE_ = grid_spec.get("Ni");
//    if (grid_spec.has("hash"))    hash_ = (std::string)grid_spec.get("hash");
//    grid_spec.get_bounding_box(bbox_);
//    grid_spec.get_points(points_);
//}

RegularLatLon::RegularLatLon( const eckit::Params& p )
{
	nsIncrement_ = p["grid_ns"];
	weIncrement_ = p["grid_ew"];

	if( ! p.get("area_s").isNil() )
	{
		bbox_ = Grid::makeBBox(p);
	}
	else
	{
		bbox_ = Grid::makeGlobalBBox();
	}

	double east  = bbox_.top_right().lon();
	double west  = bbox_.bottom_left().lon();
	double north = bbox_.top_right().lat();
	double south = bbox_.bottom_left().lat();

	nptsNS_ = computeRows( north, south, west, east );
	nptsWE_ = computeCols( west, east );

	DEBUG_VAR( nptsNS_ );
	DEBUG_VAR( nptsWE_ );

	/// @todo must solve how we do hashes, independently of the GRIB hash

	hash_ = "d0ccb07e4b36a8911817cc07539cf859"; // regular_ll 1/1

	DEBUG_HERE;
}

RegularLatLon::RegularLatLon(size_t ni, size_t nj, const Grid::BoundBox& bbox) :
	nptsNS_(nj),
	nptsWE_(ni)
{
	double east  = bbox_.top_right().lon();
	double west  = bbox_.bottom_left().lon();
	double north = bbox_.top_right().lat();
	double south = bbox_.bottom_left().lat();

	nsIncrement_ = computeIncLat();
	weIncrement_ = computeIncLon();
}

RegularLatLon::~RegularLatLon()
{
}

string RegularLatLon::hash() const
{
	 return hash_;
}

Grid::BoundBox RegularLatLon::boundingBox() const
{
	return bbox_;
}

size_t RegularLatLon::nPoints() const
{
	return nptsNS_ * nptsWE_;
}

Grid::Point RegularLatLon::latLon(size_t the_i, size_t the_j) const
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

double RegularLatLon::computeIncLat() const
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

double RegularLatLon::computeIncLon() const
{
    double east  = bbox_.top_right().lon();
    double west  = bbox_.bottom_left().lon();

    if (cols() > (east - west))
        return ((east - west)/cols());

    // Avoid truncation errors
    long inc_lon = ((east - west)/cols() + 0.5 );
    return inc_lon;
}

long RegularLatLon::computeRows(double north, double south, double west, double east) const
{
    if (north > 0.0 && south > 0.0 )
        return (north - south)/incLat() + 1;
    else
        if ( north < 0.0 && south < 0.0)
            return (fabs(north) - fabs(south))/incLat() + 1;

    return (fabs(north) + fabs(south))/incLat() + 1;
}

long RegularLatLon::computeCols(double west, double east) const
{
    return fabs((east - west)/ incLon()) + 1;
}

void RegularLatLon::coordinates( std::vector<double>& pts ) const
{
	ASSERT( pts.size() && pts.size()%2 == 0 );
	ASSERT( pts.size() == nPoints()*2 );

	const double plon = bbox_.west();    // west
	const double plat = bbox_.north();   // north;

	for( size_t j = 0; j < nptsNS_; ++j )
	{
		const double lat = plat - nsIncrement_ * j;

		for( size_t i = 0; i < nptsWE_; ++i )
		{
			const size_t idx = j*nptsWE_ + i;
			const double lon = plon + weIncrement_ * i;

			pts[ 2*idx   ] = lat;
			pts[ 2*idx+1 ] = lon;
		}
	}
}

void RegularLatLon::coordinates( std::vector<Grid::Point>& pts ) const
{
	ASSERT( pts.size() == nPoints() );

	const double plon = bbox_.west();    // west
	const double plat = bbox_.north();   // north;

	for( size_t j = 0; j < nptsNS_; ++j )
	{
		const double lat = plat - nsIncrement_ * j;

		for( size_t i = 0; i < nptsWE_; ++i )
		{
			const size_t idx = j*nptsWE_ + i;
			const double lon = plon + weIncrement_ * i;

			pts[ idx ].assign( lat, lon );
		}
	}
}

string RegularLatLon::gridType() const
{
	 return std::string("regular_ll");
}

GridSpec* RegularLatLon::spec() const
{
    GridSpec* grid_spec = new GridSpec(gridType());

    std::stringstream ss; ss << "LL" << nptsNS_ << "_" << nptsWE_;
    grid_spec->set_short_name(ss.str());

    grid_spec->set("Nj",eckit::Value(nptsNS_));
    grid_spec->set("Ni",eckit::Value(nptsWE_));

    grid_spec->set("hash",eckit::Value(hash_));
    grid_spec->set_bounding_box(bbox_);

    return grid_spec;
}

bool RegularLatLon::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const RegularLatLon&>(grid).nptsNS_ != nptsNS_) return false;
   if ( static_cast<const RegularLatLon&>(grid).nptsWE_ != nptsWE_) return false;
   if ( static_cast<const RegularLatLon&>(grid).hash_ != hash_) return false;
   if ( static_cast<const RegularLatLon&>(grid).bbox_ != bbox_) return false;

   return true;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
