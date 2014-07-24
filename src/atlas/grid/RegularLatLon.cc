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
#include "eckit/types/FloatCompare.h"
#include "eckit/config/Resource.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/RegularLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,RegularLatLon> RegularLatLon_builder("regular_ll");

RegularLatLon::RegularLatLon( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = makeBBox(p);

    DEBUG_VAR( bbox_ );
    DEBUG_VAR( bbox_.area() );

    bool built_ninj = false;
    bool built_incs = false;

    if( p.has("Nj") && p.has("Ni") )
    {
        computeGridNiNj( p["Nj"], p["Ni"] );
        built_ninj = true;
    }

    if( p.has("grid_ns_inc") && p.has("grid_ew_inc") )
    {
        computeGridIncs( p["grid_ns_inc"], p["grid_ew_inc"] );
        built_incs = true;
    }

    if( ! built_ninj && ! built_incs )
        throw BadParameter("Not enought information to build RegularLatLon", Here() );

    /// checks

	ASSERT( nptsNS_ > 1 ); // can't have a grid with just one row
	ASSERT( nptsWE_ > 1 ); // can't have a grid with just one col

    RealCompare<double> cmp( degrees_eps() );

    if( built_ninj && p.has("grid_ns_inc") )
    {
        double jInc = p["grid_ns_inc"];
        if( ! cmp( incNS_, jInc ) )
            Log::warning() << "Increment in latitude " <<  jInc << " does not match expected value " << incNS_ << std::endl;
    }

    if( built_ninj &&  p.has("grid_ew_inc") )
    {
        double iInc = p["grid_ew_inc"];
        if( ! cmp( incWE_, iInc ) )
            Log::warning() << "Increment in longitude " <<  iInc << " does not match expected value " << incWE_ << std::endl;
    }

    if( built_incs &&  p.has("Nj") )
    {
        long nj = p["Nj"];
        if( nptsNS_ != nj )
            Log::warning() << "Number of j columns " << nj << " does not match expected value " << nptsNS_ << std::endl;
    }

    if( built_incs &&  p.has("Ni") )
    {
        long ni = p["Ni"];
        if( nptsWE_ != ni )
            Log::warning() << "Number of i rows " << ni << " does not match expected value " << nptsWE_ << std::endl;
    }

    DEBUG_VAR( *spec() );
}

RegularLatLon::RegularLatLon(size_t ni, size_t nj, const Grid::BoundBox& bbox) :
	nptsNS_(nj),
	nptsWE_(ni),
	bbox_(bbox)
{
	ASSERT( nptsNS_ > 1 ); // can't have a grid with just one row
	ASSERT( nptsWE_ > 1 ); // can't have a grid with just one col

	incNS_ = computeIncLat();
	incWE_ = computeIncLon();
}

RegularLatLon::~RegularLatLon()
{
}

string RegularLatLon::uid() const
{
	std::stringstream ss;
	ss << RegularLatLon::gridTypeStr() << "_" << nptsNS_ << "_" << nptsWE_;
	return ss.str();
}

Grid::BoundBox RegularLatLon::boundingBox() const
{
	return bbox_;
}

size_t RegularLatLon::nPoints() const
{
	return nptsNS_ * nptsWE_;
}

double RegularLatLon::computeIncLat() const
{
	return (bbox_.north() - bbox_.south()) / (rows() - 1);

//	double north = bbox_.top_right().lat();
//	double south = bbox_.bottom_left().lat();

//	double north_diff_south = 0.0;
//	if (north > 0.0 && south > 0.0 ) north_diff_south = north - south;
//	else if ( north < 0.0 && south < 0.0) north_diff_south = fabs(north) - fabs(south);
//	else north_diff_south  = fabs(north) + fabs(south);

//	if (rows() > north_diff_south)
//		return north_diff_south/rows();

//	// Avoid truncation errors
//	long inc_lat = north_diff_south/(rows() + 1) + 0.5;
//	return inc_lat;
}

double RegularLatLon::computeIncLon() const
{
	return (bbox_.east() - bbox_.west()) / (cols() - 1);

//	double east  = bbox_.top_right().lon();
//	double west  = bbox_.bottom_left().lon();

//	if (cols() > (east - west))
//		return ((east - west)/cols());

//	// Avoid truncation errors
//	long inc_lon = ((east - west)/cols() + 0.5 );
    //	return inc_lon;
}

long RegularLatLon::computeRows() const
{
    return (bbox_.east() - bbox_.west()) / incLon() + 1;
}

long RegularLatLon::computeCols() const
{
    return (bbox_.north() - bbox_.south()) / incLat() + 1;
}

void RegularLatLon::computeGridNiNj(long Ni, long Nj)
{
    nptsNS_ = Ni;
    nptsWE_ = Nj;

    incNS_ = computeIncLat();
    incWE_ = computeIncLon();
}

void RegularLatLon::computeGridIncs(double incNS, double incWE)
{
    incNS_ = incNS;
    incWE_ = incWE;

    nptsNS_ = computeRows();
    nptsWE_ = computeCols();
}

//long RegularLatLon::computeRows(double north, double south, double west, double east) const
//{
//    if (north > 0.0 && south > 0.0 )
//        return (north - south)/incLat() + 1;
//    else
//        if ( north < 0.0 && south < 0.0)
//            return (fabs(north) - fabs(south))/incLat() + 1;

//    return (fabs(north) + fabs(south))/incLat() + 1;
//}

//long RegularLatLon::computeCols(double west, double east) const
//{
//    return fabs((east - west)/ incLon()) + 1;
//}

void RegularLatLon::coordinates( std::vector<double>& pts ) const
{
	ASSERT( pts.size() && pts.size()%2 == 0 );
	ASSERT( pts.size() == nPoints()*2 );

	const double plon = bbox_.west();    // west
	const double plat = bbox_.north();   // north;

	for( size_t j = 0; j < nptsNS_; ++j )
	{
		const double lat = plat - incNS_ * j;

		for( size_t i = 0; i < nptsWE_; ++i )
		{
			const size_t idx = j*nptsWE_ + i;
			const double lon = plon + incWE_ * i;

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
		const double lat = plat - incNS_ * j;

		for( size_t i = 0; i < nptsWE_; ++i )
		{
			const size_t idx = j*nptsWE_ + i;
			const double lon = plon + incWE_ * i;

			pts[ idx ].assign( lat, lon );
		}
	}
}

string RegularLatLon::gridType() const
{
	return RegularLatLon::gridTypeStr();
}

GridSpec* RegularLatLon::spec() const
{
	GridSpec* grid_spec = new GridSpec( gridType() );

	grid_spec->uid( uid() );

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
