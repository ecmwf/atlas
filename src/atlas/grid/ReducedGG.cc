/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "grib_api.h"

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

ConcreteBuilderT1<Grid,ReducedGG> ReducedGG_builder("reduced_gg");

ReducedGG::ReducedGG( const eckit::Params& p )
{
	if( !p.get("hash").isNil() )
		hash_ = p["hash"].as<std::string>();

	bbox_ = makeBBox(p);

	gaussN_ = p["GaussN"];
	nj_     = p["Nj"];

	ASSERT( nj_ == 2*gaussN_ );

	ValueList nlats = p["NLats"];
	rgSpec_.resize(nlats.size());
	for( size_t i = 0; i < nlats.size(); ++i)
		rgSpec_[i] = nlats[i];

	ASSERT( rgSpec_.size() == 2 * gaussN_ ); // number of lines of latitude should, be twice the gaussN_

	latitudes_.resize( rgSpec_.size() );

	/// @todo this code should be moved into Atlas library and co-maintained with NA section
	grib_get_gaussian_latitudes(gaussN_, &latitudes_[0]);

	/// @todo this should not be necessary here -- check it...
	///       we should either use latitudes (angle from equator) or colatitutes (angle from pole)
	///       and class ReducedGG should stick to that definition
	//	if (jScansPositively_ == 1 )
	//	   std::reverse(latitudes_.begin(), latitudes_.end());

	computePoints();

	nbDataPoints_ = p["nbDataPoints"];

	DEBUG_VAR( nbDataPoints_ );
	DEBUG_VAR( points_.size() );

	ASSERT( points_.size() == nbDataPoints_ );
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
	for( size_t i = 0; i < points_.size(); ++i )
	{
		pts[ 2*i   ] = points_[i].lat();
		pts[ 2*i+1 ] = points_[i].lon();
	}
}

void ReducedGG::coordinates(std::vector<Grid::Point>& pts) const
{
	ASSERT( pts.size() == nPoints() );
	std::copy(points_.begin(),points_.end(),pts.begin());
}

string ReducedGG::gridType() const
{
	return ReducedGG::gridTypeStr();
}

GridSpec* ReducedGG::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType());

   grid_spec->uid( uid() );
   grid_spec->set("gaussianNumber",eckit::Value(gaussN_));

   grid_spec->set("hash",eckit::Value(hash_));

   grid_spec->set_bounding_box(bbox_);
   grid_spec->set_rgspec(rgSpec_);
   grid_spec->set_latitudes(latitudes_);

   return grid_spec;
}

bool ReducedGG::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const ReducedGG&>(grid).gaussN_ != gaussN_) return false;
   if ( static_cast<const ReducedGG&>(grid).hash_ != hash_) return false;
   if ( static_cast<const ReducedGG&>(grid).bbox_ != bbox_) return false;
   if ( static_cast<const ReducedGG&>(grid).rgSpec_ != rgSpec_) return false;
   if ( static_cast<const ReducedGG&>(grid).latitudes_ != latitudes_) return false;
   if ( static_cast<const ReducedGG&>(grid).points_ != points_) return false;

   return true;
}

void ReducedGG::computePoints()
{
	ASSERT( points_.size() == 0 );

	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west = bbox_.west();
	const double east = bbox_.east();

//	DEBUG_VAR(north);
//	DEBUG_VAR(south);
//	DEBUG_VAR(west);
//	DEBUG_VAR(east);

	const std::vector<double>& lats = latitudes_;

	for ( size_t i = 0;  i < rgSpec_.size(); ++i )
	{
		// check latitudes bound box
		if( ( lats[i] <= north && lats[i] >= south )
			|| isEqual(lats[i],north)
			|| isEqual(lats[i],south) )
		{
			const long nlats = rgSpec_[i];

			ASSERT( nlats > 0 );

			const double delta_lat = 360.0/nlats;
			double plon = 0;

			for( long k = 0; k < nlats; ++k )
			{
				if( ( plon >= west && plon <= east ) || isEqual(plon,west) || isEqual(plon,east) )
				{
					points_.push_back( Grid::Point( lats[i], plon ) );
				}
//				else { std::cout << "skipped " << Grid::Point( lats[i], plon ) << std::endl; }

				plon += delta_lat;
			}
		}
//		else { std::cout << "skipped lats[" << i << "] : " << lats[i] << std::endl; }
	}
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
