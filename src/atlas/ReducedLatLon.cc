/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/geometry/Point2.h"
#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/value/Value.h"

#include "atlas/GridSpec.h"
#include "atlas/ReducedLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {


//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,ReducedLatLon> ReducedLatLon_builder( ReducedLatLon::gridTypeStr() );

ReducedLatLon::ReducedLatLon( const eckit::Params& p ) :
	npts_(0),
	nsInc_(0),
	nptsNS_(0)
{
   if( !p.get("hash").isNil() )
      hash_ = p["hash"].as<std::string>();

   bbox_ = makeBBox(p);

   if( p.has("Nj") )
   {
      nptsNS_ = p["Nj"];
   }

   if( p.has("grid_lat_inc") )
   {
	  nsInc_ = p["grid_lat_inc"];
   }

   if( p.has("NPtsPerLat") )
   {
      ValueList nlats = p["NPtsPerLat"];
      nbPtsPerLat_.resize(nlats.size());
      for( size_t i = 0; i < nlats.size(); ++i)
         nbPtsPerLat_[i] = nlats[i];
   }
   else
   {
      computeNPtsPerLat(nbPtsPerLat_);
   }

   npts_ = computeNPts();
}

ReducedLatLon::~ReducedLatLon()
{
}

string ReducedLatLon::uid() const
{
	std::stringstream ss;
	ss << ReducedLatLon::gridTypeStr() << "_" << nptsNS_;
	return ss.str();
}

struct ReducedLatLon_CoordDD
{
	ReducedLatLon_CoordDD( std::vector<double>& pts ) : pts_(pts) {}
	void operator()(double lat, double lon)
	{
		pts_.push_back( lat );
		pts_.push_back( lon );
	}
	std::vector<double>& pts_;
};


void ReducedLatLon::coordinates(std::vector<double>& r ) const
{
	r.clear();
	r.reserve( nPoints() * 2);

	ReducedLatLon_CoordDD f(r);

	iterate(f);
}

struct ReducedLatLon_Coord
{
	ReducedLatLon_Coord( std::vector<Grid::Point>& pts ) : pts_(pts) {}
	void operator()(double lat, double lon)
	{
		pts_.push_back( Grid::Point(lat,lon) );
	}
	std::vector<Grid::Point>& pts_;
};


void ReducedLatLon::coordinates( std::vector<Grid::Point>& pts) const
{
	pts.clear();

	pts.reserve( nPoints() );

	ReducedLatLon_Coord f(pts);

	iterate(f);
}

size_t ReducedLatLon::nPoints() const
{
	return npts_;
}

string ReducedLatLon::gridType() const
{
	return ReducedLatLon::gridTypeStr();
}

GridSpec ReducedLatLon::spec() const
{
   GridSpec grid_spec(gridType());

   grid_spec.uid( uid() );

   grid_spec.set("Nj",eckit::Value(nptsNS_));
   grid_spec.set("grid_lat_inc",eckit::Value(nsInc_));

   grid_spec.set("hash",eckit::Value(hash_));

   grid_spec.set_bounding_box(bbox_);
   grid_spec.set_npts_per_lat(nbPtsPerLat_);

   return grid_spec;
}

bool ReducedLatLon::same(const Grid& grid) const
{
	return spec() == grid.spec();
}

long ReducedLatLon::computeRows() const
{
	return (bbox_.north() - bbox_.south()) / incLat() + 1;
}

struct ReducedLatLon_Counter
{
	ReducedLatLon_Counter() : count_(0) {}
	void operator()(double lat, double lon)
	{
		++count_;
	}
	size_t count_;
};

size_t ReducedLatLon::computeNPts() const
{
	ReducedLatLon_Counter f;

	iterate(f);

	return f.count_;
}

void ReducedLatLon::computeNPtsPerLat( std::vector<long>& )
{
   // Not clear how this is computed.
   // Note: For guassain number, we used the pre-defined tabulated values.
   NOTIMP;
}

template <class T>
void ReducedLatLon::iterate( T& f ) const
{
	RealCompare<double> isEqual(degrees_eps());

	const double north = bbox_.north();
	const double south = bbox_.south();
	const double west = bbox_.west();
	const double east = bbox_.east();

	double plat = north;
	for( size_t j = 0; j < nbPtsPerLat_.size() ; ++j) {

	   long no_of_points_along_latitude = nbPtsPerLat_[j];
	   if (no_of_points_along_latitude > 0 ) {

		  if ( isEqual(plat,north) || isEqual(plat,south) || ( plat < north && plat > south)) {

			 double east_west_grid_length = 360.0/no_of_points_along_latitude;
			 double plon = 0;
			 for(int k = 0; k < no_of_points_along_latitude; k++) {

				if ( (isEqual(plon,west) || isEqual(plon,east)) || ( plon < east && plon > west)) {

				   ASSERT(plat <= 90.0 && plat >= -90.0);
				   ASSERT(plon < 360.0 && plon >= 0);
				   f(plat,plon);
				}
				plon += east_west_grid_length;
				geometry::reduceTo2Pi(plon);

			 }
		  }
	   }
	   plat -= nsInc_;
	}
}

//-----------------------------------------------------------------------------


} // namespace eckit
