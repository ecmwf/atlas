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
#include "eckit/types/FloatCompare.h"
#include "eckit/value/Value.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/RegularGG.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,RegularGG> RegularGG_builder( RegularGG::gridTypeStr() );

RegularGG::RegularGG( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = makeBBox(p);

	gaussN_ = p["GaussN"];

	ASSERT( gaussN_ > 1 );

	if( p.has("Nj") )
	{
		long Nj = p["Nj"];
		if( Nj != nj() )
			Log::warning() << "Number of points along meridian" << Nj << " does not match expected value " << nj() << std::endl;
	}

	ni_ = p["Ni"];

	ASSERT( ni_ > 1 );

	nbDataPoints_ = p["nbDataPoints"];
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

void RegularGG::coordinates( std::vector<double>& pts ) const
{
	ASSERT( pts.size() && pts.size()%2 == 0 );
	ASSERT( pts.size() == nPoints()*2 );

	std::vector<double> latitudes( nj() );
	grib_get_gaussian_latitudes(gaussN_, &latitudes[0]);

	std::vector<Grid::Point> points;
	computePoints(latitudes,points);

	for( size_t i = 0; i < points.size(); ++i )
	{
		pts[ 2*i   ] = points[i].lat();
		pts[ 2*i+1 ] = points[i].lon();
	}
}

void RegularGG::coordinates(std::vector<Grid::Point>& pts) const
{
	ASSERT( pts.size() == nbDataPoints_ );

	std::vector<double> latitudes( nj() );
	grib_get_gaussian_latitudes(gaussN_, &latitudes[0]);

	computePoints(latitudes,pts);
}

string RegularGG::gridType() const
{
	return RegularGG::gridTypeStr();
}

GridSpec RegularGG::spec() const
{
   GridSpec grid_spec(gridType());

   grid_spec.uid( uid() );

   grid_spec.set("GaussN", gaussN_);

   grid_spec.set("Ni", ni_);
   grid_spec.set("Nj", nj());

   grid_spec.set("grid_lon_inc", computeIncLon() );

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

	pts.resize( nPoints() );

	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west  = bbox_.west();
	const double east  = bbox_.east();

	const double delta_lon = computeIncLon();

	size_t n = 0;
	for ( size_t j = 0;  j < nj(); ++j )
	{
		// check latitudes bound box
		if( ( lats[j] <= north && lats[j] >= south ) || isEqual(lats[j],north) || isEqual(lats[j],south) )
		{
			double plon = 0;

			for( long k = 0; k < ni_; ++k )
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

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
