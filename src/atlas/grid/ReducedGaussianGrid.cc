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
#include "atlas/grid/ReducedGaussianGrid.h"
#include "atlas/grid/GribRead.h"

using namespace eckit;
using namespace std;


namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ? No.
// NPoints: is this just the grids points, or includes area points, of area does not fit grid
//          assumes it is grid points inclusive of the area.

ReducedGaussianGrid::ReducedGaussianGrid(grib_handle* handle)
: GribGrid(handle),
  gaussianNumber_(0),
  nj_(0)
{
   Log::info() << "Build a ReducedGaussianGrid  " << std::endl;

   // Extract the guassian grid attributes from the grib handle

   GRIB_CHECK(grib_get_long(handle,"numberOfParallelsBetweenAPoleAndTheEquator",&gaussianNumber_),0);
   GRIB_CHECK(grib_get_long(handle,"Nj",&nj_),0);

   // get reduced grid specification. These are number of point's along the lines of latitude
   size_t rgSpecLength = 0;
   GRIB_CHECK(grib_get_size(handle,"pl",&rgSpecLength),0);
   rgSpec_.resize(rgSpecLength);
   GRIB_CHECK(grib_get_long_array(handle,"pl",&rgSpec_[0],&rgSpecLength),0);


   // Need to check AREA geometry, which uses scanning mode ???
   // .......

   // This provides the 'y' co-ordinate of each line of latitude
   latitudes_.resize(2 *gaussianNumber_ );
   grib_get_gaussian_latitudes(gaussianNumber_, &latitudes_[0]);

   // number of lines of latitude should, be twice the numberOfParallelsBetweenAPoleAndTheEquator
   ASSERT( rgSpec_.size() == 2 *gaussianNumber_);


   // Create point list based on area. To avoid rounding errors determine if we have
   // Global area, then use simple algorithm.
   if (isGlobalNorthSouth() && isGlobalWestEast()) {
      Log::info() << " GLOBAL " << std::endl;
      for ( int i = 0; i < 2*gaussianNumber_ && i < rgSpec_.size(); i++ ) {
         long no_of_points_along_latitude = rgSpec_[i];
         double east_west_grid_length = 360.0/no_of_points_along_latitude;
         double plon = 0;
         for(int k = 0; k < no_of_points_along_latitude; k++) {
            points_.push_back( Point( latitudes_[i], plon ) );
            plon += east_west_grid_length;
         }
      }
   }
   else {
      Log::info() << " LOCAL " << std::endl;
      for ( int i = 0; i < 2*gaussianNumber_ && i < rgSpec_.size(); i++ ) {
         if ( FloatCompare::is_equal(latitudes_[i],north_,epsilon()) ) {
            // Log::info() << " same north latitudes_[i] " <<  latitudes_[i] << " north_ " << north_ << std::endl;
            add_point(i);
            continue;
         }
         if ( FloatCompare::is_equal(latitudes_[i],south_,epsilon()) ) {
            add_point(i);
            continue;
         }
         if ( latitudes_[i] < north_ && latitudes_[i] > south_) {
            add_point(i);
         }
      }
   }
   if (jScansPositively_ == 1 )
      std::reverse(latitudes_.begin(), latitudes_.end());

   Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
   Log::info() << " epsilon()                                      " << epsilon() << std::endl;
   Log::info() << " gaussianNumber_                                " << gaussianNumber_ << std::endl;
   Log::info() << " nj                                             " << nj_ << std::endl;
   Log::info() << " isGlobalNorthSouth()                           " << isGlobalNorthSouth() << std::endl;
   Log::info() << " isGlobalWestEast()                             " << isGlobalWestEast() << std::endl;
   Log::info() << " pl-size                                        " << rgSpecLength <<  " pl[0]="<< rgSpec_[0] << " pl[1]=" << rgSpec_[1] << std::endl;
   Log::info() << " latitudes_.size ()                             " << latitudes_.size() << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively_ << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively_ << std::endl;
   Log::info() << " scanning_mode                                  " << Grid::scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << east_ << std::endl;
   Log::info() << " EXPECTED longitudeOfLastGridPointInDegrees     " << 360.0 - (360.0/(gaussianNumber_*4.0)) << std::endl;
   Log::info() << " numberOfDataPoints                             " << numberOfDataPoints_ << std::endl;
   int no_of_points_in_pl = 0;
   for(int i = 0; i < rgSpec_.size(); i++) no_of_points_in_pl += rgSpec_[i];
   Log::info() << " no_of_points_in_pl                             " << no_of_points_in_pl << std::endl;
   Log::info() << " points_.size()                                 " << points_.size() << std::endl;

   ASSERT(points_.size() == numberOfDataPoints_);

   // Check point list compared with grib
   GribRead::comparePointList(points_,epsilon(),handle);
}

ReducedGaussianGrid::~ReducedGaussianGrid()
{
    Log::info() << "Destroy a ReducedGaussianGrid" << std::endl;
}

void ReducedGaussianGrid::add_point(int lat_index)
{
   long no_of_points_along_latitude = rgSpec_[lat_index];
   double east_west_grid_length = 360.0/no_of_points_along_latitude;
   double plon = 0;
   for(int k = 0; k < no_of_points_along_latitude; k++) {
      if ( (plon > west_ && plon < east_) ||
           (FloatCompare::is_equal(plon,west_,epsilon())  ||
            FloatCompare::is_equal(plon,east_,epsilon()) )
          ) {
         points_.push_back( Point( latitudes_[lat_index], plon ) );
         plon += east_west_grid_length;
      }
   }
}

Grid::BoundBox ReducedGaussianGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Point(south_,west_),Point(north_,east_));
}

void ReducedGaussianGrid::coordinates( Grid::Coords& r ) const
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

bool ReducedGaussianGrid::isGlobalNorthSouth() const
{
   return (gaussianNumber_*2 == nj_);
   //std::cout << "isGlobalNorthSouth " << (fabs(south_ - latitudes_[gaussianNumber_*2-1])) << "     " << fabs( north_ - latitudes_[0]) << "\n";
   //return (fabs(south_ - latitudes_[gaussianNumber_*2-1])) <= AREA_FACTOR && fabs( north_ - latitudes_[0]) <= AREA_FACTOR;
}

bool ReducedGaussianGrid::isGlobalWestEast() const
{
   /// AREA_FACTOR is added because grib has precision for 3 dec places.
   double res = east_ - west_ + 90.0 / gaussianNumber_ + AREA_FACTOR;
//   cout << " ReducedGaussianGrid::isGlobalWestEast() 90.0 / gaussianNumber_ = " << double(90.0 /(double)gaussianNumber_) << endl;
//   cout << " ReducedGaussianGrid::isGlobalWestEast() double(360.0 /((double)4*gaussianNumber_)  = " << (360.0/((double)4*gaussianNumber_)) << endl;
//   cout << " ReducedGaussianGrid::isGlobalWestEast() RES = " << res << endl;
   return res > 360.0 || FloatCompare::is_equal(res,360.0,epsilon());
}
//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
