/*
 * (C) Copyright 1996-2014 ECMWF.
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

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/ReducedLatLon.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

ConcreteBuilderT1<Grid,ReducedLatLon> ReducedLatLon_builder( ReducedLatLon::gridTypeStr() );

ReducedLatLon::ReducedLatLon( const eckit::Params& p )
 : nsIncrement_(0),nptsNS_(0)
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
      nsIncrement_ = p["grid_lat_inc"];
   }

   if( p.has("NPtsPerLat") )
   {
      ValueList nlats = p["NPtsPerLat"];
      nbPtsPerLat_.resize(nlats.size());
      for( size_t i = 0; i < nlats.size(); ++i)
         nbPtsPerLat_[i] = nlats[i];
   }
   else {
      computeNPtsPerLat(nbPtsPerLat_);
   }
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

void ReducedLatLon::coordinates(std::vector<double>& r ) const
{
   if (r.capacity() == 0)  r.reserve( nPoints() * 2);
   r.clear();

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

                  ASSERT(plat < 90.0 && plat > -90.0);
                  ASSERT(plon < 360.0 && plon >= 0);
                  r.push_back( plat );
                  r.push_back( plon );
               }
               plon += east_west_grid_length;
            }
         }
      }
      plat -= nsIncrement_;
   }
}

void ReducedLatLon::coordinates(std::vector<Grid::Point>& pts) const
{
   if (pts.size() == 0) {
      pts.resize( nPoints() );
   }

   RealCompare<double> isEqual(degrees_eps());

   const double north = bbox_.north();
   const double south = bbox_.south();
   const double west = bbox_.west();
   const double east = bbox_.east();

   int i = 0;
   double plat = north;
   for( size_t j = 0; j < nbPtsPerLat_.size() ; ++j) {

      long no_of_points_along_latitude = nbPtsPerLat_[j];
      if (no_of_points_along_latitude > 0 ) {

         if ( isEqual(plat,north) || isEqual(plat,south) || ( plat < north && plat > south)) {

            double east_west_grid_length = 360.0/no_of_points_along_latitude;
            double plon = 0;
            for(int k = 0; k < no_of_points_along_latitude; k++) {

               if ( (isEqual(plon,west) || isEqual(plon,east)) || ( plon < east && plon > west)) {

                  ASSERT(plat < 90.0 && plat > -90.0);
                  ASSERT(plon < 360.0 && plon >= 0);
                  ASSERT(i < pts.size());
                  pts[i].assign( plat, plon );
                  i++;
               }
               plon += east_west_grid_length;
            }
         }
      }
      plat -= nsIncrement_;
   }
}

size_t ReducedLatLon::nPoints() const
{
   size_t nbDataPoints = 0;

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

                  ASSERT(plat < 90.0 && plat > -90.0);
                  ASSERT(plon < 360.0 && plon >= 0);
                  nbDataPoints++;
               }
               plon += east_west_grid_length;
            }
         }
      }
      plat -= nsIncrement_;
   }
//
//   for( size_t j = 0; j < nbPtsPerLat_.size(); ++j) {
//      nbDataPoints += nbPtsPerLat_[j];
//   }
   return nbDataPoints;
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
   grid_spec.set("grid_lat_inc",eckit::Value(nsIncrement_));

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

void ReducedLatLon::computeNPtsPerLat( std::vector<long>& )
{
   // Not clear how this is computed.
   // Note: For guassain number, we used the pre-defined tabulated values.
   NOTIMP;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
