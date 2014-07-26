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

ConcreteBuilderT1<Grid,RegularLatLon> RegularLatLon_builder( RegularLatLon::gridTypeStr() );

RegularLatLon::RegularLatLon( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = makeBBox(p);

	bool built_ninj = false;
    bool built_incs = false;

	if( p.has("Nj") && p.has("Ni") )
	{
		computeGridNiNj( p["Nj"], p["Ni"] );
		built_ninj = true;
	}

	if( !built_ninj && p.has("grid_lat_inc") && p.has("grid_lon_inc") )
    {
		computeGridIncs( p["grid_lat_inc"], p["grid_lon_inc"] );
        built_incs = true;
    }

    if( ! built_ninj && ! built_incs )
        throw BadParameter("Not enought information to build RegularLatLon", Here() );

    /// checks

	ASSERT( nptsNS_ > 1 ); // can't have a grid with just one row
	ASSERT( nptsWE_ > 1 ); // can't have a grid with just one col

    RealCompare<double> cmp( degrees_eps() );

	if( built_ninj && p.has("grid_lat_inc") )
    {
		double jInc = p["grid_lat_inc"];
        if( ! cmp( incNS_, jInc ) )
            Log::warning() << "Increment in latitude " <<  jInc << " does not match expected value " << incNS_ << std::endl;
    }

	if( built_ninj &&  p.has("grid_lon_inc") )
    {
		double iInc = p["grid_lon_inc"];
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

	if( p.has("nbDataPoints") )
	{
		long npts = p["nbDataPoints"];
		if( npts != nPoints() )
			Log::warning() << "Number of data points " << npts << " does not match expected value " << nPoints() << std::endl;
	}

//	DEBUG_VAR( bbox_ );
//  DEBUG_VAR( bbox_.area() );
//	DEBUG_VAR( built_ninj );
//	DEBUG_VAR( built_incs );
//	DEBUG_VAR( spec() );
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
}

double RegularLatLon::computeIncLon() const
{
	return (bbox_.east() - bbox_.west()) / (cols() - 1);
}

long RegularLatLon::computeRows() const
{
	return (bbox_.north() - bbox_.south()) / incLat() + 1;
}

long RegularLatLon::computeCols() const
{
	return (bbox_.east() - bbox_.west()) / incLon() + 1;
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

GridSpec RegularLatLon::spec() const
{
	GridSpec grid_spec( gridType() );

	grid_spec.uid( uid() );

	grid_spec.set("Nj", nptsNS_ );
	grid_spec.set("Ni", nptsWE_ );

	grid_spec.set("grid_lat_inc", incNS_);
	grid_spec.set("grid_lon_inc", incWE_);

	grid_spec.set("hash", hash_ );

	grid_spec.set_bounding_box(bbox_);

    return grid_spec;
}

bool RegularLatLon::same(const Grid& grid) const
{
	return spec() == grid.spec();
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
