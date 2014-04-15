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
#include "atlas/grid/RegularGaussianGrid.h"

using namespace eckit;

const double epsilon = 1.25e-10;
bool zero(const double v) { return abs(v) < epsilon; }
bool same(const double a, const double b) { return zero(a-b); }
bool zero(const double v, double epsilon) { return abs(v) < epsilon; }
bool same(const double a, const double b, double epsilon) { return zero(a-b, epsilon); }

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

RegularGaussianGrid::RegularGaussianGrid(grib_handle* handle)
: gaussianNumber_(0),
  north_(0.0),south_(0.0),west_(0.0),east_(0.0)
{
   Log::info() << "Build a RegularGaussianGrid  " << std::endl;

   // Extract the guassian grid attributes from the grib handle
   if (handle == NULL) throw std::runtime_error("NULL grib_handle");
   hash_ = grib_hash(handle);

   long iScansNegatively = 0, jScansPositively = 0;
   GRIB_CHECK(grib_get_long(handle,"iScansNegatively",&iScansNegatively),0);
   GRIB_CHECK(grib_get_long(handle,"jScansPositively",&jScansPositively),0);
   int scanning_mode = Grid::scanningMode(iScansNegatively,jScansPositively);

   GRIB_CHECK(grib_get_double(handle,"latitudeOfFirstGridPointInDegrees",&north_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfFirstGridPointInDegrees",&west_),0);
   GRIB_CHECK(grib_get_double(handle,"latitudeOfLastGridPointInDegrees",&south_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfLastGridPointInDegrees",&east_),0);

   GRIB_CHECK(grib_get_long(handle,"numberOfParallelsBetweenAPoleAndTheEquator",&gaussianNumber_),0);


//   GRIB_CHECK(grib_get_double(handle,"jDirectionIncrementInDegrees",&nsIncrement_),0);
//   GRIB_CHECK(grib_get_double(handle,"iDirectionIncrementInDegrees",&weIncrement_),0);
//
//   GRIB_CHECK(grib_get_long(handle,"Nj",&nptsNS_),0);
//   GRIB_CHECK(grib_get_long(handle,"Ni",&nptsWE_),0);

   long nb_nodes = 0;
   grib_get_long(handle,"numberOfDataPoints",&nb_nodes);

   // Need to check AREA geometry, which uses scanning mode ???
   // .......

   double array[2 *gaussianNumber_];
   grib_get_gaussian_latitudes(gaussianNumber_, array);


   for ( int i = 0; i < 2*gaussianNumber_; i++ )
   {
      if ( same(array[i],north_, 10e-2) ) {
         latitudes_.push_back(array[i]);
         continue;
      }
      if ( same(array[i],south_, 10e-2) ) {
         latitudes_.push_back(array[i]);
         continue;
      }
      if ( array[i] < north_ && array[i] > south_)
         latitudes_.push_back(array[i]);
   }

//   if (jScansPositively == 1 )
//      std::reverse(latitudes_.begin(), latitudes_.end());

   // These needs to be reviewed. TODO
   long nptsWE = 4 * gaussianNumber_ ;
   long weIncrement = 90./gaussianNumber_ ; // ??
   double plon = west_;
   double plat = north_;
   for(size_t i = 0 ; i < latitudes_.size(); i++) {

      for( size_t j = 0; j < nptsWE; ++j) {
         points_.push_back( Point( latitudes_[i], plon ) );
         plon += weIncrement;
      }
   }

   long nptsNS = 2 * gaussianNumber_ ;

   Log::info() << " gaussianNumber_                                " << gaussianNumber_ << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively << std::endl;
   Log::info() << " scanning_mode                                  " << scanning_mode << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << east_ << std::endl;
   Log::info() << " iDirectionIncrementInDegrees(west-east   incr) " << weIncrement << std::endl;
   Log::info() << " nptsNS                                         " << nptsNS << std::endl;
   Log::info() << " nptsWE                                         " << nptsWE << std::endl;
   Log::info() << " numberOfDataPoints                             " << nb_nodes << std::endl;
   Log::info() << " -----------------------------------------------" << std::endl;
   Log::info() << " points_.size()  " << points_.size() << "       nb_nodes " << nb_nodes << std::endl << std::endl;

   ASSERT(points_.size() == nb_nodes);
}

RegularGaussianGrid::~RegularGaussianGrid()
{
    Log::info() << "Destroy a RegularGaussianGrid" << std::endl;
}

Grid::Point RegularGaussianGrid::latLon(size_t the_i, size_t the_j) const
{
    long nptsWE = 4 * gaussianNumber_ ;
    long weIncrement = 90./gaussianNumber_ ; // ??
    double plon = west_;
    for(size_t i = 0 ; i < latitudes_.size(); i++) {

       for( size_t j = 0; j < nptsWE; ++j) {
          if ( i== the_i && j== the_j) {
             return  Point( latitudes_[i], plon );
          }
          plon += weIncrement;
       }
    }

   return Grid::Point();
}

std::string RegularGaussianGrid::hash() const
{
    return hash_;
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


//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
