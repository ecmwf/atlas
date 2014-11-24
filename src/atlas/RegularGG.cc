/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/atlas_config.h"

#ifdef ECKIT_HAVE_GRIB
#include "grib_api.h" // for grib_get_gaussian_latitudes()
#endif

#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/value/Value.h"

#include "atlas/GridSpec.h"
#include "atlas/RegularGG.h"

using namespace eckit;
using namespace std;

namespace atlas {


//------------------------------------------------------------------------------------------------------

RegularGG::RegularGG( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = make_bounding_box(p);

	gaussN_ = p["N"];

	ASSERT( gaussN_ > 1 );

	if( p.has("Nj") )
	{
		long Nj = p["Nj"];
		if( Nj != nj() )
			Log::warning() << "Number of points along meridian" << Nj << " does not match expected value " << nj() << std::endl;
	}

	if( p.has("Ni") )
	{
		ni_ = p["Ni"];
		ASSERT( ni_ == 4*gaussN_ );
	}
	else
	{
		ni_ = 4*gaussN_;
		Log::warning() << "Assuming number of points along parallel to be 4 * GaussianNumber " << ni_ << std::endl;
	}

	if( p.has("npts") )
		npts_ = p["npts"];
	else
	{
		std::vector<double> latitudes;
		computeLatitudes(latitudes);
		npts_ = computeNPoints(latitudes);
	}

	ASSERT( npts_ > 0 );
}

RegularGG::~RegularGG()
{
}

string RegularGG::uid() const
{
	std::stringstream ss;
	ss << gridTypeStr() << "_" << gaussN_;
	return ss.str();
}

void RegularGG::lonlat( double pts[] ) const
{
	std::vector<double> latitudes;
	computeLatitudes(latitudes);

	std::vector<Grid::Point> points;
	computePoints(latitudes,points);

  int c(0);
	for( size_t i = 0; i < points.size(); ++i )
	{
		pts[c++] = points[i].lon();
		pts[c++] = points[i].lat();
	}
}

void RegularGG::lonlat(std::vector<Grid::Point>& pts) const
{
	ASSERT( pts.size() == npts_ );

	std::vector<double> latitudes;
	computeLatitudes(latitudes);

	computePoints(latitudes,pts);
}

string RegularGG::grid_type() const
{
	return RegularGG::gridTypeStr();
}

GridSpec RegularGG::spec() const
{
   GridSpec grid_spec(grid_type());

   grid_spec.uid( uid() );

   grid_spec.set("N", gaussN_);

   grid_spec.set("Ni", ni_);
   grid_spec.set("Nj", nj());

   grid_spec.set("lon_inc", computeIncLon() );

   grid_spec.set("hash",hash_);
   grid_spec.set_bounding_box(bbox_);

   return grid_spec;
}

bool RegularGG::same(const Grid& grid) const
{
	return spec() == grid.spec();
}

double RegularGG::computeIncLon() const
{
	return 360.0/ni_;
}

void RegularGG::computePoints( const std::vector<double>& lats, std::vector<Point>& pts ) const
{
	ASSERT( lats.size() == nj() );

	pts.resize( npts() );

	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west  = bbox_.west();
	const double east  = bbox_.east();

	const double delta_lon = computeIncLon();

	size_t n = 0;
	for ( size_t j = 0;  j < nj(); ++j )
	{
		if( ( lats[j] <= north && lats[j] >= south ) || isEqual(lats[j],north) || isEqual(lats[j],south) )
		{
			double plon = 0;

			for( long k = 0; k < ni_; ++k )
			{
				if( ( plon >= west && plon <= east ) || isEqual(plon,west) || isEqual(plon,east) )
				{
					pts[n].assign( plon, lats[j] );
					++n;
				}
				plon += delta_lon;
			}
		}
	}

	ASSERT( n == npts_ );
}

long RegularGG::computeNPoints(const std::vector<double>& lats) const
{
	ASSERT( lats.size() == nj() );

	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west  = bbox_.west();
	const double east  = bbox_.east();

	const double delta_lon = computeIncLon();

	size_t n = 0;
	for ( size_t j = 0;  j < nj(); ++j )
	{
		if( ( lats[j] <= north && lats[j] >= south ) || isEqual(lats[j],north) || isEqual(lats[j],south) )
		{
			double plon = 0;

			for( long k = 0; k < ni_; ++k )
			{
				if( ( plon >= west && plon <= east ) || isEqual(plon,west) || isEqual(plon,east) )
				{
					++n;
				}
				plon += delta_lon;
			}
		}
	}
	return n;
}

void RegularGG::computeLatitudes(std::vector<double>& lats) const
{
	lats.resize( nj() );

	/// @todo this should not be necessary here -- check it...
	///       we should either use latitudes (angle from equator) or colatitutes (angle from pole)
	///       and class ReducedGG should stick to that definition
	//	if (jScansPositively_ == 1 )
	//	   std::reverse(lats.begin(), lats.end());

#ifdef ECKIT_HAVE_GRIB
	/// @todo this code should be moved into Atlas library and co-maintained with NA section
	grib_get_gaussian_latitudes(gaussN_, &lats[0]);
#else
	NOTIMP;
#endif
}

//-----------------------------------------------------------------------------


} // namespace eckit
