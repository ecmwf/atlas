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
#include "eckit/types/FloatCompare.h"
#include "eckit/grib/GribAccessor.h"
#include "atlas/grid/RegularGaussianGrid.h"
#include "atlas/grid/GribRead.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

RegularGaussianGrid::RegularGaussianGrid(grib_handle* handle)
: GribGrid(handle),
  gaussianNumber_(0),
  nj_(0)
{
   Log::info() << "Build a RegularGaussianGrid  " << std::endl;

   // Extract the guassian grid attributes from the grib handle

   GRIB_CHECK(grib_get_long(handle,"numberOfParallelsBetweenAPoleAndTheEquator",&gaussianNumber_),0);
   GRIB_CHECK(grib_get_long(handle,"Nj",&nj_),0);


   // Need to check AREA geometry, which uses scanning mode ???
   // .......

   double array[2 *gaussianNumber_];
   grib_get_gaussian_latitudes(gaussianNumber_, array);


   double nptsWE = 4 * gaussianNumber_ ;
   double weIncrement = 360.0/nptsWE;
   if (isGlobalNorthSouth() && isGlobalWestEast()) {
      Log::info() << " GLOBAL                      weIncrement = " << weIncrement << std::endl;
      for ( int i = 0; i < 2*gaussianNumber_; i++ ) {
         double lat = array[i];
         ASSERT(lat < 90.0 && lat > -90.0);
         latitudes_.push_back(array[i]);
      }
      if (jScansPositively_ == 1 )
         std::reverse(latitudes_.begin(), latitudes_.end());
      for(size_t i = 0 ; i < latitudes_.size(); i++) {
         double plon = 0;
         for( size_t j = 0; j < nptsWE; ++j) {
            ASSERT(latitudes_[i] < 90.0);
            points_.push_back( Point( latitudes_[i], plon ) );
            plon += weIncrement;
         }
      }
   }
   else {
      Log::info() << " LOCAL " << std::endl;
      for ( int i = 0; i < 2*gaussianNumber_; i++ ) {
         if ( FloatCompare::is_equal(array[i],north_,epsilon()) ) {
            latitudes_.push_back(array[i]);
            continue;
         }
         if ( FloatCompare::is_equal(array[i],south_,epsilon()) ) {
            latitudes_.push_back(array[i]);
            continue;
         }
         if ( array[i] < north_ && array[i] > south_)
            latitudes_.push_back(array[i]);
      }
      if (jScansPositively_ == 1 )
          std::reverse(latitudes_.begin(), latitudes_.end());
      for(size_t i = 0 ; i < latitudes_.size(); i++) {
         double plon = west_;
         for( size_t j = 0; j < nptsWE; ++j) {
            if (plon < east_ || FloatCompare::is_equal(plon,east_)) {
               points_.push_back( Point( latitudes_[i], plon ) );
               plon += weIncrement;
            }
         }
      }
   }

   Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
   Log::info() << " epsilon()                                      " << epsilon() << std::endl;
   Log::info() << " gaussianNumber_                                " << gaussianNumber_ << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively_ << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively_ << std::endl;
   Log::info() << " scanning_mode                                  " << Grid::scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << east_ << std::endl;
   Log::info() << " iDirectionIncrementInDegrees(west-east   incr) " << weIncrement << std::endl;
   Log::info() << " nptsNS                                         " << 2 * gaussianNumber_ << std::endl;
   Log::info() << " nptsWE                                         " << nptsWE << std::endl;
   Log::info() << " numberOfDataPoints                             " << numberOfDataPoints_ << std::endl;
   Log::info() << " -----------------------------------------------" << std::endl;
   Log::info() << " points_.size()  " << points_.size() << "       numberOfDataPoints " << numberOfDataPoints_ << std::endl;
   Log::info() << " point[0]                               " << points_[0].lat() << ", " << points_[0].lon() <<  std::endl;
   Log::info() << " point[1]                               " << points_[1].lat() << ", " << points_[1].lon() <<  std::endl;

   ASSERT(points_.size() == numberOfDataPoints_);

   // Check point list compared with grib
   GribRead::comparePointList(points_,epsilon(),handle);
}

RegularGaussianGrid::~RegularGaussianGrid()
{
    Log::info() << "Destroy a RegularGaussianGrid" << std::endl;
}

Grid::Point RegularGaussianGrid::latLon(size_t the_i, size_t the_j) const
{
    long nptsWE = 4 * gaussianNumber_ ;
    long weIncrement = 360.0 / nptsWE;
    for(size_t i = 0 ; i < latitudes_.size(); i++) {

       double plon = west_;
       for( size_t j = 0; j < nptsWE; ++j) {
          if ( i== the_i && j== the_j) {
             return  Point( latitudes_[i], plon );
          }
          plon += weIncrement;
       }
    }

   return Grid::Point();
}

Grid::BoundBox RegularGaussianGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Point(south_,west_),Point(north_,east_));
}

void RegularGaussianGrid::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

/// AREA_FACTOR is added because GRIB has precision for 3 dec places.
/// For instance east for N640 is 359.8593750 intstead of 359.859
static const double AREA_FACTOR = 1.e-3;

bool RegularGaussianGrid::isGlobalNorthSouth() const
{
   return (gaussianNumber_*2 == nj_);
   //std::cout << "isGlobalNorthSouth " << (fabs(south_ - latitudes_[gaussianNumber_*2-1])) << "     " << fabs( north_ - latitudes_[0]) << "\n";
   //return (fabs(south_ - latitudes_[gaussianNumber_*2-1])) <= AREA_FACTOR && fabs( north_ - latitudes_[0]) <= AREA_FACTOR;
}

bool RegularGaussianGrid::isGlobalWestEast() const
{
   /// AREA_FACTOR is added because grib has precision for 3 dec places.
   double res = east_ - west_ + 90.0 / gaussianNumber_ + AREA_FACTOR;
//   cout << " RegularGaussianGrid::isGlobalWestEast() 90.0 / gaussianNumber_ = " << double(90.0 /(double)gaussianNumber_) << endl;
//   cout << " RegularGaussianGrid::isGlobalWestEast() double(360.0 /((double)4*gaussianNumber_)  = " << (360.0/((double)4*gaussianNumber_)) << endl;
//   cout << " RegularGaussianGrid::isGlobalWestEast() RES = " << res << endl;
   return res > 360.0 || FloatCompare::is_equal(res,360.0,epsilon());
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
