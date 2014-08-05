/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "grib_api.h" // for grib_get_gaussian_latitudes()

#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/value/Value.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/ReducedGG.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,ReducedGG> ReducedGG_builder( ReducedGG::gridTypeStr() );

ReducedGG::ReducedGG( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = makeBBox(p);

	gaussN_ = p["GaussN"];

	ASSERT( gaussN_ > 1 );

	if( p.has("Nj") )
	{
		long nj = p["Nj"];
		if( nj != 2*gaussN_ )
			Log::warning() << "Number of j columns " << nj << " does not match expected value " << 2*gaussN_ << std::endl;
	}

	if( p.has("NPtsPerLat") )
	{
		ValueList nlats = p["NPtsPerLat"];
		rgSpec_.resize(nlats.size());
		for( size_t i = 0; i < nlats.size(); ++i)
			rgSpec_[i] = nlats[i];
	}
	else
		computeNPtsPerLat(rgSpec_);

	ASSERT( rgSpec_.size() == 2 * gaussN_ ); // number of lines of latitude should, be twice the gaussN_

	if( p.has("NPtsPerLat") )
		nbDataPoints_ = p["nbDataPoints"];
	else
	{
		std::vector<double>  latitudes;
		computeLatitues(latitudes);
		nbDataPoints_ = computeNPoints(latitudes);
	}

	ASSERT( nbDataPoints_ > 0 );
}

ReducedGG::ReducedGG(long gaussN) : gaussN_(gaussN)
{
	ASSERT( gaussN_ > 1 );

	/// @todo needs to have computeNPtsPerLat() implemented

	NOTIMP;
}

ReducedGG::~ReducedGG()
{
}

string ReducedGG::uid() const
{
	std::stringstream ss;
	ss << gridTypeStr() << "_" << gaussN_;
	return ss.str();
}

void ReducedGG::coordinates( std::vector<double>& pts ) const
{
	ASSERT( pts.size() && pts.size()%2 == 0 );
	ASSERT( pts.size() == nPoints()*2 );

	std::vector<double>  latitudes;
	computeLatitues(latitudes);

	std::vector<Grid::Point> points;
	computePoints(latitudes,points);

	for( size_t i = 0; i < points.size(); ++i )
	{
		pts[ 2*i   ] = points[i].lat();
		pts[ 2*i+1 ] = points[i].lon();
	}
}

void ReducedGG::coordinates(std::vector<Grid::Point>& pts) const
{
	ASSERT( pts.size() == nbDataPoints_ );

	std::vector<double>  latitudes;
	computeLatitues(latitudes);

	computePoints(latitudes,pts);
}

string ReducedGG::gridType() const
{
	return ReducedGG::gridTypeStr();
}

GridSpec ReducedGG::spec() const
{
   GridSpec grid_spec(gridType());

   grid_spec.uid( uid() );
   grid_spec.set("GaussN", gaussN_);

   grid_spec.set("hash", hash_);

   grid_spec.set_bounding_box(bbox_);
   grid_spec.set_rgspec(rgSpec_);

   return grid_spec;
}

bool ReducedGG::same(const Grid& grid) const
{
	return spec() == grid.spec();
}

void ReducedGG::computeLatitues(std::vector<double>& lats) const
{
	lats.resize(rgSpec_.size());

	/// @todo this should not be necessary here -- check it...
	///       we should either use latitudes (angle from equator) or colatitutes (angle from pole)
	///       and class ReducedGG should stick to that definition
	//	if (jScansPositively_ == 1 )
	//	   std::reverse(lats.begin(), lats.end());

	/// @todo this code should be moved into Atlas library and co-maintained with NA section
	grib_get_gaussian_latitudes(gaussN_, &lats[0]);
}

void ReducedGG::computePoints( const std::vector<double>& lats, std::vector<Point>& pts ) const
{
	ASSERT( lats.size() == rgSpec_.size() );

	pts.resize( nbDataPoints_ );

	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west = bbox_.west();
	const double east = bbox_.east();

	size_t n = 0;
	for ( size_t j = 0;  j < rgSpec_.size(); ++j )
	{
		if( ( lats[j] <= north && lats[j] >= south ) || isEqual(lats[j],north) || isEqual(lats[j],south) )
		{
			const long npts_per_lat = rgSpec_[j];

			ASSERT( npts_per_lat > 0 );

			const double delta_lon = 360.0/npts_per_lat;
			double plon = 0;

			for( long k = 0; k < npts_per_lat; ++k )
			{
				if( ( plon >= west && plon <= east ) || isEqual(plon,west) || isEqual(plon,east) )
				{
					pts[n].assign( lats[j], plon );
					++n;
				}
				plon += delta_lon;
			}
		}
	}

	ASSERT( n == nbDataPoints_ );
}

long ReducedGG::computeNPoints( const std::vector<double>& lats ) const
{
	ASSERT( lats.size() == rgSpec_.size() );

	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west = bbox_.west();
	const double east = bbox_.east();

	long n = 0;
	for ( size_t j = 0;  j < rgSpec_.size(); ++j )
	{
		if( ( lats[j] <= north && lats[j] >= south ) || isEqual(lats[j],north) || isEqual(lats[j],south) )
		{
			const long npts_per_lat = rgSpec_[j];

			ASSERT( npts_per_lat > 0 );

			const double delta_lon = 360.0/npts_per_lat;
			double plon = 0;

			for( long k = 0; k < npts_per_lat; ++k )
			{
				if( ( plon >= west && plon <= east ) || isEqual(plon,west) || isEqual(plon,east) )
				{
					++n;
				}
				plon += delta_lon;
			}
		}
	}
  /* nawd: no return implemented here? */ NOTIMP;
}


void ReducedGG::computeNPtsPerLat(std::vector<long>& nlats)
{
	NOTIMP;

	/// @todo deduce nlats from gauss number, probably using tablated values
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
