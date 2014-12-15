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

#include "atlas/GridSpec.h"
#include "atlas/grids/RegularLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

RegularLatLon::RegularLatLon( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = make_bounding_box(p);

	bool built_ninj = false;
    bool built_incs = false;

	if( p.has("Nj") && p.has("Ni") )
	{
		computeGridNiNj( p["Nj"], p["Ni"] );
		built_ninj = true;
	}

	if( !built_ninj && p.has("lat_inc") && p.has("lon_inc") )
    {
		computeGridIncs( p["lat_inc"], p["lon_inc"] );
        built_incs = true;
    }

    if( ! built_ninj && ! built_incs )
        throw BadParameter("Not enought information to build RegularLatLon", Here() );

    /// checks

	ASSERT( nptsNS_ > 1 ); // can't have a grid with just one row
	ASSERT( nptsWE_ > 1 ); // can't have a grid with just one col

    RealCompare<double> cmp( degrees_eps() );

	if( built_ninj && p.has("lat_inc") )
    {
		double jInc = p["lat_inc"];
        if( ! cmp( incNS_, jInc ) )
            Log::warning() << "Increment in latitude " <<  jInc << " does not match expected value " << incNS_ << std::endl;
    }

	if( built_ninj &&  p.has("lon_inc") )
    {
		double iInc = p["lon_inc"];
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

	if( p.has("npts") )
	{
		long nb_pts = p["npts"];
		if( nb_pts != npts() )
			Log::warning() << "Number of data points " << nb_pts << " does not match expected value " << npts() << std::endl;
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
	ss << RegularLatLon::grid_type_str() << "_" << nptsNS_ << "_" << nptsWE_;
	return ss.str();
}

Grid::BoundBox RegularLatLon::bounding_box() const
{
	return bbox_;
}

size_t RegularLatLon::npts() const
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

void RegularLatLon::lonlat( double pts[] ) const
{
	const double plon = bbox_.west();    // west
	const double plat = bbox_.north();   // north;

	for( size_t j = 0; j < nptsNS_; ++j )
	{
		const double lat = plat - incNS_ * j;

		for( size_t i = 0; i < nptsWE_; ++i )
		{
			const size_t idx = j*nptsWE_ + i;
			const double lon = plon + incWE_ * i;

			pts[ 2*idx   ] = lon;
			pts[ 2*idx+1 ] = lat;
		}
	}
}

void RegularLatLon::lonlat( std::vector<Grid::Point>& pts ) const
{
  pts.resize(npts());
  const double plon = bbox_.west();    // west
	const double plat = bbox_.north();   // north;

	for( size_t j = 0; j < nptsNS_; ++j )
	{
		const double lat = plat - incNS_ * j;

		for( size_t i = 0; i < nptsWE_; ++i )
		{
			const size_t idx = j*nptsWE_ + i;
			const double lon = plon + incWE_ * i;

			pts[ idx ].assign( lon, lat );
		}
	}
}

string RegularLatLon::grid_type() const
{
	return RegularLatLon::grid_type_str();
}

GridSpec RegularLatLon::spec() const
{
	GridSpec grid_spec( grid_type() );

	grid_spec.uid( uid() );

	grid_spec.set("Nj", nptsNS_ );
	grid_spec.set("Ni", nptsWE_ );

	grid_spec.set("lat_inc", incNS_);
	grid_spec.set("lon_inc", incWE_);

	grid_spec.set("hash", hash_ );

	grid_spec.set_bounding_box(bbox_);

    return grid_spec;
}

bool RegularLatLon::same(const Grid& grid) const
{
	return spec() == grid.spec();
}

//-----------------------------------------------------------------------------


} // namespace grids
} // namespace atlas
