/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

// ==================================================================================
// gribs use the following convention: (from Shahram)
//
// Horizontally:  Points scan in the +i (+x) direction
// Vertically:    Points scan in the -j (-y) direction
//
// The way I verified this was to look at our SAMPLE files (which IFS uses).
// I also verified that IFS does not modify the scanning modes
// so whatever the samples say, is the convention
// ==================================================================================

#include <stdexcept>

#include "eckit/log/Log.h"
#include "eckit/grib/GribAccessor.h"
#include "atlas/grid/GribRead.h"
#include "atlas/grid/RegularLatLonGrid.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

RegularLatLonGrid::RegularLatLonGrid(grib_handle* handle)
: nsIncrement_(0),weIncrement_(0),nptsNS_(0),nptsWE_(0),
  north_(0.0),south_(0.0),west_(0.0),east_(0.0),editionNumber_(0)
{
   Log::info() << "Build a RegularLatLonGrid  " << std::endl;

   // Extract the regluar lat long grid attributes from the grib handle
   if (handle == NULL) throw std::runtime_error("NULL grib_handle");

   GRIB_CHECK(grib_get_long(handle,"editionNumber",&editionNumber_),0);

   hash_ = grib_hash(handle);

   long iScansNegatively = 0, jScansPositively = 0;
   GRIB_CHECK(grib_get_long(handle,"iScansNegatively",&iScansNegatively),0);
   GRIB_CHECK(grib_get_long(handle,"jScansPositively",&jScansPositively),0);

   GRIB_CHECK(grib_get_double(handle,"latitudeOfFirstGridPointInDegrees",&north_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfFirstGridPointInDegrees",&west_),0);
   GRIB_CHECK(grib_get_double(handle,"latitudeOfLastGridPointInDegrees",&south_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfLastGridPointInDegrees",&east_),0);

   GRIB_CHECK(grib_get_double(handle,"jDirectionIncrementInDegrees",&nsIncrement_),0);
   GRIB_CHECK(grib_get_double(handle,"iDirectionIncrementInDegrees",&weIncrement_),0);

   GRIB_CHECK(grib_get_long(handle,"Nj",&nptsNS_),0);
   GRIB_CHECK(grib_get_long(handle,"Ni",&nptsWE_),0);

   long nb_nodes = 0;
   grib_get_long(handle,"numberOfDataPoints",&nb_nodes);

   int scanning_mode = Grid::scanningMode(iScansNegatively,jScansPositively);


   // Need to check AREA geometry, which uses scanning mode ???
   // .......

   double plat = north_;
   points_.reserve( (nptsNS_ + 1) * (nptsWE_ + 1) );
   for( size_t j = 0; j < nptsNS_; ++j) {
      double plon = west_;
      for( size_t i = 0; i < nptsWE_; ++i) {
         points_.push_back( Point( plat, plon ) );
         plon += weIncrement_;
      }
      plat -= nsIncrement_;
   }

   Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively << std::endl;
   Log::info() << " scanning_mode                                  " << scanning_mode << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << east_ << std::endl;
   Log::info() << " jDirectionIncrementInDegrees(north-south incr) " << nsIncrement_ << std::endl;
   Log::info() << " iDirectionIncrementInDegrees(west-east   incr) " << weIncrement_ << std::endl;
   Log::info() << " Nj(num of points North South)                  " << nptsNS_ << std::endl;
   Log::info() << " Ni(num of points West East)                    " << nptsWE_ << std::endl;
   Log::info() << " numberOfDataPoints                             " << nb_nodes << std::endl;
   Log::info() << " -----------------------------------------------" << std::endl;
   Log::info() << " computeIncLat() " << computeIncLat() << "      nsIncrement_ " << nsIncrement_ << std::endl;
   Log::info() << " computeIncLon() " << computeIncLon() << "      weIncrement_ " << nsIncrement_ << std::endl;
   Log::info() << " computeRows()   " << computeRows(north_,south_,west_,east_) << "     nptsNS_ " << nptsNS_ << std::endl;
   Log::info() << " computeCols()   " << computeCols(west_,east_) <<  "     nptsWE_ " << nptsWE_ << std::endl;
   Log::info() << " points_.size()  " << points_.size() << "       nb_nodes " << nb_nodes << std::endl << std::endl;

   ASSERT(nsIncrement_ == computeIncLat());
   ASSERT(weIncrement_ == computeIncLon());
   ASSERT(nptsNS_ == computeRows(north_,south_,west_,east_));
   ASSERT(nptsWE_ == computeCols(west_,east_));
   ASSERT(points_.size() == nb_nodes);

   // Check point list compared with grib
   GribRead::comparePointList(points_,epsilon(),handle);
}

RegularLatLonGrid::~RegularLatLonGrid()
{
    Log::info() << "Destroy a RegularLatLonGrid" << std::endl;
}

Grid::Point RegularLatLonGrid::latLon(size_t the_i, size_t the_j) const
{
   double plon = west_;
   double plat = north_;
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

std::string RegularLatLonGrid::hash() const
{
    return hash_;
}

Grid::BoundBox RegularLatLonGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Point(south_,west_),Point(north_,east_));
}

void RegularLatLonGrid::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

long RegularLatLonGrid::computeIncLat() const
{
   double north_diff_south = 0.0;
   if (north_ > 0.0 && south_ > 0.0 ) north_diff_south = north_ - south_;
   else if ( north_ < 0.0 && south_ < 0.0) north_diff_south = fabs(north_) - fabs(south_);
   else north_diff_south  = fabs(north_) + fabs(south_);

   return ( north_diff_south/(rows() + 1) + 0.5);
}

long RegularLatLonGrid::computeIncLon() const
{
   return ((east_ - west_)/cols() + 0.5 );
}

long RegularLatLonGrid::computeRows(double north, double south, double west, double east) const
{
   if (north > 0.0 && south > 0.0 ) return (north - south)/nsIncrement_ + 1;
   else if ( north < 0.0 && south < 0.0) return (fabs(north) - fabs(south))/nsIncrement_ + 1;

   return (fabs(north) + fabs(south))/nsIncrement_ + 1;
}

long RegularLatLonGrid::computeCols(double west, double east) const
{
   return fabs((east - west)/weIncrement_) + 1;
}

double RegularLatLonGrid::epsilon() const
{
   // grib edition 1 - milli-degrees
   // grib edition 2 - micro-degrees or could be defined by the keys: "subdivisionsOfBasicAngle" and "basicAngleOfTheInitialProductionDomain"
   // Therefore the client needs access to this when dealing with double based comparisons (for tolerances)
   return (editionNumber_ == 1) ? 1e-3 : 1e-6;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
