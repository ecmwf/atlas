/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include <stdexcept>

#include "atlas/grid/GridBuilder.h"
#include "atlas/grid/ReducedGaussianGrid.h"
#include "atlas/grid/RegularGaussianGrid.h"
#include "atlas/grid/RegularLatLonGrid.h"
#include "atlas/grid/Unstructured.h"

#include "eckit/log/Log.h"
#include "eckit/geometry/Point3.h"
#include "eckit/grib/GribAccessor.h"
#include "eckit/types/FloatCompare.h"

using namespace std;
using namespace eckit;
using namespace eckit::geometry;


/// AREA_FACTOR is added because GRIB has precision for 3 dec places.
/// For instance east for N640 is 359.8593750 intstead of 359.859
static const double AREA_FACTOR = 1.e-3;

namespace atlas {
namespace grid {

typedef std::vector< atlas::grid::Grid::Point > PointList;
static PointList* read_number_of_data_points(grib_handle *h);


//=====================================================================================
GridBuilder::GridBuilder(){}
GridBuilder::~GridBuilder(){}

//=====================================================================================

GribGridBuilder::GribGridBuilder() {}

GribGridBuilder::~GribGridBuilder() {}

Grid::Ptr GribGridBuilder::build(eckit::LocalPathName pathname) const
{
   FILE* fp = fopen(pathname.c_str(),"r");
   if (!fp) {
      std::stringstream ss; ss << "Could not open file " << pathname.c_str();
      throw std::runtime_error(ss.str());
   }

   int err;
   grib_handle* handle = grib_handle_new_from_file(0,fp,&err);
   if (err != 0 || handle == 0) {
      std::stringstream ss; ss << "Could not create grib handle for file " << pathname.c_str() << "  err=" << err;
      throw std::runtime_error(ss.str());
   }

   return build_grid_from_grib_handle(handle);
}


Grid::Ptr GribGridBuilder::build_grid_from_grib_handle( grib_handle* handle ) const
{
   char string_value[64];
   size_t len = sizeof(string_value)/sizeof(char);
   if (grib_get_string(handle,"gridType",string_value,&len)  != 0) {
      throw std::runtime_error("grib_get_string failed for gridType") ;
   }

   if (strncasecmp(string_value,"regular_ll",10) == 0) {
      // Custom arguments for Lat long grid
      GribRegularLatLonGrid maker( handle );
      return Grid::Ptr(new RegularLatLonGrid(maker.hash(),
                                             maker.boundingBox(),
                                             maker.coordinates(),
                                             maker.incLat(),  // nsIncrement
                                             maker.incLon(),  // weIncrement
                                             maker.rows(),    // nptsNS
                                             maker.cols()     // nptsWE
      ));
   }
   //   else if (strncasecmp(string_value,"sh",2) == 0) {
   //      GribSphericalHarmonicGrid maker( handle );
   //      return SphericalHarmonicGrid( maker.hash(), maker.boundingBox(),...);
   //   }
   //   else if (strncasecmp(string_value,"reduced_ll",10) == 0) {
   //      GribReducedLatLonGrid maker( handle );
   //      return ReducedLatLonGrid( maker.hash(), maker.boundingBox(),...);
   //   }
   else if (strncasecmp(string_value,"reduced_gg",10) == 0) {
      GribReducedGaussianGrid maker( handle );
      return Grid::Ptr(new ReducedGaussianGrid(maker.hash(),
                                               maker.boundingBox(),
                                               maker.coordinates(),
                                               maker.latitudes(),
                                               maker.gaussianNumber()));
   }
   else if (strncasecmp(string_value,"regular_gg",10) == 0) {
      GribRegularGaussianGrid maker( handle );
      return Grid::Ptr(new RegularGaussianGrid(maker.hash(),
                                               maker.boundingBox(),
                                               maker.coordinates(),
                                               maker.latitudes(),
                                               maker.gaussianNumber()));
   }

   // Unknown grid type, get extract data points form the grib handle
   return Grid::Ptr(new grid::Unstructured( read_number_of_data_points(handle), grib_hash(handle) ));
}


GribGridBuilder& GribGridBuilder::instance()
{
   static GribGridBuilder* obj = 0;

   if( !obj )
      obj = new GribGridBuilder();

   return *obj;
}

void GribGridBuilder::known_grid_types(std::set<std::string>& grids)
{
   grids.insert("regular_ll");
   grids.insert("reduced_ll");
   grids.insert("mercator");
   grids.insert("lambert");
   grids.insert("polar_stereographic");
   grids.insert("UTM");
   grids.insert("simple_polyconic");
   grids.insert("albers");
   grids.insert("miller");
   grids.insert("rotated_ll");
   grids.insert("stretched_ll");
   grids.insert("stretched_rotated_ll");
   grids.insert("regular_gg");
   grids.insert("rotated_gg");
   grids.insert("stretched_gg");
   grids.insert("stretched_rotated_gg");
   grids.insert("reduced_gg");
   grids.insert("sh");
   grids.insert("rotated_sh");
   grids.insert("stretched_sh");
   grids.insert("stretched_rotated_sh");
   grids.insert("space_view");
   grids.insert("unknown");
   grids.insert("unknown_PLPresent");
}

// ==================================================================================================

GribGridMaker::GribGridMaker(grib_handle* handle)
: editionNumber_(0),
  north_(0.0),south_(0.0),west_(0.0),east_(0.0),
  epsilon_(1e-6),
  numberOfDataPoints_(0)
{
   if (handle == NULL)
      throw std::runtime_error("NULL grib_handle");

   GRIB_CHECK(grib_get_long(handle,"editionNumber",&editionNumber_),0);
   epsilon_ = (editionNumber_ == 1) ? 1e-3 : 1e-6;

   hash_ = grib_hash(handle);

   GRIB_CHECK(grib_get_long(handle,"iScansNegatively",&iScansNegatively_),0);
   GRIB_CHECK(grib_get_long(handle,"jScansPositively",&jScansPositively_),0);

   GRIB_CHECK(grib_get_double(handle,"latitudeOfFirstGridPointInDegrees",&north_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfFirstGridPointInDegrees",&west_),0);
   GRIB_CHECK(grib_get_double(handle,"latitudeOfLastGridPointInDegrees",&south_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfLastGridPointInDegrees",&east_),0);

   // Check area
   ASSERT(north_ > south_);
   ASSERT(north_ < 90.0  || FloatCompare::is_equal(north_,90.0,epsilon_));
   ASSERT(south_ < 90.0  || FloatCompare::is_equal(south_,90.0,epsilon_));
   ASSERT(north_ > -90.0 || FloatCompare::is_equal(north_,-90.0,epsilon_));
   ASSERT(south_ > -90.0 || FloatCompare::is_equal(south_,-90.0,epsilon_));

   // make sure we are in range 0-360.0 ??
   while(west_ < 0) west_ = west_ + 360.0;
   while(east_ < 0) east_ = east_ + 360.0;
   while(west_ > 360.0) west_ = west_ - 360.0;
   while(east_ > 360.0) east_ = east_ - 360.0;

   grib_get_long(handle,"numberOfDataPoints",&numberOfDataPoints_);
}

GribGridMaker::~GribGridMaker(){}


int GribGridMaker::scanningMode(long iScansNegatively, long jScansPositively)
{
   if(iScansNegatively == 0 && jScansPositively == 0)
      return 1;
   if(iScansNegatively == 0 && jScansPositively == 1)
      return 2;
   if(iScansNegatively == 1 && jScansPositively == 0)
      return 3;
   if(iScansNegatively == 1 && jScansPositively == 1)
      return 4;

   ASSERT(iScansNegatively == 0 || iScansNegatively == 1);
   ASSERT(jScansPositively == 0 || jScansPositively == 1);

   return 1;
}

void GribGridMaker::comparePointList(const std::vector<Grid::Point>& points,  double epsilon, grib_handle* handle)
{
   // Check point list compared with grib
   PointList* pointlist1 = read_number_of_data_points(handle);
   const std::vector<atlas::grid::Grid::Point>& pntlist = *pointlist1;
   ASSERT( points.size() == pntlist.size());
   int print_point_list = 0;
   std::streamsize old_precision = cout.precision();
   for(size_t i =0;  i < points.size() && i < pntlist.size(); i++) {
      if (!FloatCompare::is_equal(points[i].lat(),pntlist[i].lat(),epsilon) ||
               !FloatCompare::is_equal(points[i].lon(),pntlist[i].lon(),epsilon))
      {
         if (print_point_list == 0) Log::info() << " Point list DIFFER, show first 10, epsilon(" << epsilon << ")" << std::endl;
         Log::info() << setw(3) << i << " :"
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << points[i].lat() << ", "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << points[i].lon() << "  "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << pntlist[i].lat() << ", "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << pntlist[i].lon()
                  << std::endl;
         print_point_list++;
         if (print_point_list > 10) break;
      }
   }

   // reset precision
   Log::info() << std::setprecision(old_precision);

   delete pointlist1;
}

// ================================================================================================

GribReducedGaussianGrid::GribReducedGaussianGrid(grib_handle* handle)
: GribGridMaker(handle),
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
            points_.push_back( Grid::Point( latitudes_[i], plon ) );
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
   Log::info() << " scanning_mode                                  " << scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
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
   comparePointList(points_,epsilon(),handle);
}

GribReducedGaussianGrid::~GribReducedGaussianGrid()
{
   Log::info() << "Destroy a ReducedGaussianGrid" << std::endl;
}

void GribReducedGaussianGrid::add_point(int lat_index)
{
   long no_of_points_along_latitude = rgSpec_[lat_index];
   double east_west_grid_length = 360.0/no_of_points_along_latitude;
   double plon = 0;
   for(int k = 0; k < no_of_points_along_latitude; k++) {
      if ( (plon > west_ && plon < east_) ||
               (FloatCompare::is_equal(plon,west_,epsilon())  ||
                        FloatCompare::is_equal(plon,east_,epsilon()) )
      ) {
         points_.push_back( Grid::Point( latitudes_[lat_index], plon ) );
         plon += east_west_grid_length;
      }
   }
}

Grid::BoundBox GribReducedGaussianGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox( Grid::Point(south_,west_),Grid::Point(north_,east_));
}

void GribReducedGaussianGrid::coordinates( Grid::Coords& r ) const
{
   ASSERT( r.size() == points_.size() );

   for( size_t i = 0; i < points_.size(); ++i )
   {
      r.lat(i) = points_[i].lat();
      r.lon(i) = points_[i].lon();
   }
}

bool GribReducedGaussianGrid::isGlobalNorthSouth() const
{
   return (gaussianNumber_*2 == nj_);
   //std::cout << "isGlobalNorthSouth " << (fabs(south_ - latitudes_[gaussianNumber_*2-1])) << "     " << fabs( north_ - latitudes_[0]) << "\n";
   //return (fabs(south_ - latitudes_[gaussianNumber_*2-1])) <= AREA_FACTOR && fabs( north_ - latitudes_[0]) <= AREA_FACTOR;
}

bool GribReducedGaussianGrid::isGlobalWestEast() const
{
   /// AREA_FACTOR is added because grib has precision for 3 dec places.
   double res = east_ - west_ + 90.0 / gaussianNumber_ + AREA_FACTOR;
   //   cout << " ReducedGaussianGrid::isGlobalWestEast() 90.0 / gaussianNumber_ = " << double(90.0 /(double)gaussianNumber_) << endl;
   //   cout << " ReducedGaussianGrid::isGlobalWestEast() double(360.0 /((double)4*gaussianNumber_)  = " << (360.0/((double)4*gaussianNumber_)) << endl;
   //   cout << " ReducedGaussianGrid::isGlobalWestEast() RES = " << res << endl;
   return res > 360.0 || FloatCompare::is_equal(res,360.0,epsilon());
}

// ========================================================================================

GribRegularGaussianGrid::GribRegularGaussianGrid(grib_handle* handle)
: GribGridMaker(handle),
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
            points_.push_back( Grid::Point( latitudes_[i], plon ) );
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
               points_.push_back( Grid::Point( latitudes_[i], plon ) );
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
   Log::info() << " scanning_mode                                  " << scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
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
   comparePointList(points_,epsilon(),handle);
}

GribRegularGaussianGrid::~GribRegularGaussianGrid()
{
   Log::info() << "Destroy a RegularGaussianGrid" << std::endl;
}

Grid::Point GribRegularGaussianGrid::latLon(size_t the_i, size_t the_j) const
{
   long nptsWE = 4 * gaussianNumber_ ;
   long weIncrement = 360.0 / nptsWE;
   for(size_t i = 0 ; i < latitudes_.size(); i++) {

      double plon = west_;
      for( size_t j = 0; j < nptsWE; ++j) {
         if ( i== the_i && j== the_j) {
            return  Grid::Point( latitudes_[i], plon );
         }
         plon += weIncrement;
      }
   }

   return Grid::Point();
}

Grid::BoundBox GribRegularGaussianGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Grid::Point(south_,west_),Grid::Point(north_,east_));
}

void GribRegularGaussianGrid::coordinates( Grid::Coords& r ) const
{
   ASSERT( r.size() == points_.size() );

   for( size_t i = 0; i < points_.size(); ++i )
   {
      r.lat(i) = points_[i].lat();
      r.lon(i) = points_[i].lon();
   }
}

bool GribRegularGaussianGrid::isGlobalNorthSouth() const
{
   return (gaussianNumber_*2 == nj_);
   //std::cout << "isGlobalNorthSouth " << (fabs(south_ - latitudes_[gaussianNumber_*2-1])) << "     " << fabs( north_ - latitudes_[0]) << "\n";
   //return (fabs(south_ - latitudes_[gaussianNumber_*2-1])) <= AREA_FACTOR && fabs( north_ - latitudes_[0]) <= AREA_FACTOR;
}

bool GribRegularGaussianGrid::isGlobalWestEast() const
{
   /// AREA_FACTOR is added because grib has precision for 3 dec places.
   double res = east_ - west_ + 90.0 / gaussianNumber_ + AREA_FACTOR;
   //   cout << " RegularGaussianGrid::isGlobalWestEast() 90.0 / gaussianNumber_ = " << double(90.0 /(double)gaussianNumber_) << endl;
   //   cout << " RegularGaussianGrid::isGlobalWestEast() double(360.0 /((double)4*gaussianNumber_)  = " << (360.0/((double)4*gaussianNumber_)) << endl;
   //   cout << " RegularGaussianGrid::isGlobalWestEast() RES = " << res << endl;
   return res > 360.0 || FloatCompare::is_equal(res,360.0,epsilon());
}

// =====================================================================================

GribRegularLatLonGrid::GribRegularLatLonGrid(grib_handle* handle)
: GribGridMaker(handle),
  nsIncrement_(0),weIncrement_(0),nptsNS_(0),nptsWE_(0)
{
   Log::info() << "Build a RegularLatLonGrid  " << std::endl;

   // Extract the regluar lat long grid attributes from the grib handle

   GRIB_CHECK(grib_get_double(handle,"jDirectionIncrementInDegrees",&nsIncrement_),0);
   GRIB_CHECK(grib_get_double(handle,"iDirectionIncrementInDegrees",&weIncrement_),0);

   GRIB_CHECK(grib_get_long(handle,"Nj",&nptsNS_),0);
   GRIB_CHECK(grib_get_long(handle,"Ni",&nptsWE_),0);

   double plat = north_;
   points_.reserve( (nptsNS_ + 1) * (nptsWE_ + 1) );
   for( size_t j = 0; j < nptsNS_; ++j) {
      double plon = west_;
      for( size_t i = 0; i < nptsWE_; ++i) {
         points_.push_back( Grid::Point( plat, plon ) );
         plon += weIncrement_;
      }
      plat -= nsIncrement_;
   }

   Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively_ << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively_ << std::endl;
   Log::info() << " scanning_mode                                  " << scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << east_ << std::endl;
   Log::info() << " jDirectionIncrementInDegrees(north-south incr) " << nsIncrement_ << std::endl;
   Log::info() << " iDirectionIncrementInDegrees(west-east   incr) " << weIncrement_ << std::endl;
   Log::info() << " Nj(num of points North South)                  " << nptsNS_ << std::endl;
   Log::info() << " Ni(num of points West East)                    " << nptsWE_ << std::endl;
   Log::info() << " numberOfDataPoints                             " << numberOfDataPoints_ << std::endl;
   Log::info() << " -----------------------------------------------" << std::endl;
   Log::info() << " computeIncLat() " << computeIncLat() << "      nsIncrement_ " << nsIncrement_ << std::endl;
   Log::info() << " computeIncLon() " << computeIncLon() << "      weIncrement_ " << nsIncrement_ << std::endl;
   Log::info() << " computeRows()   " << computeRows(north_,south_,west_,east_) << "     nptsNS_ " << nptsNS_ << std::endl;
   Log::info() << " computeCols()   " << computeCols(west_,east_) <<  "     nptsWE_ " << nptsWE_ << std::endl;
   Log::info() << " points_.size()  " << points_.size() << "       numberOfDataPoints_ " << numberOfDataPoints_ << std::endl << std::endl;

   ASSERT(nsIncrement_ == computeIncLat());
   ASSERT(weIncrement_ == computeIncLon());
   ASSERT(nptsNS_ == computeRows(north_,south_,west_,east_));
   ASSERT(nptsWE_ == computeCols(west_,east_));
   ASSERT(points_.size() == numberOfDataPoints_);

   // Check point list compared with grib
   comparePointList(points_,epsilon(),handle);
}

GribRegularLatLonGrid::~GribRegularLatLonGrid()
{
   Log::info() << "Destroy a RegularLatLonGrid" << std::endl;
}

Grid::Point GribRegularLatLonGrid::latLon(size_t the_i, size_t the_j) const
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

Grid::BoundBox GribRegularLatLonGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Grid::Point(south_,west_),Grid::Point(north_,east_));
}

void GribRegularLatLonGrid::coordinates( Grid::Coords& r ) const
{
   ASSERT( r.size() == points_.size() );

   for( size_t i = 0; i < points_.size(); ++i )
   {
      r.lat(i) = points_[i].lat();
      r.lon(i) = points_[i].lon();
   }
}

long GribRegularLatLonGrid::computeIncLat() const
{
   double north_diff_south = 0.0;
   if (north_ > 0.0 && south_ > 0.0 ) north_diff_south = north_ - south_;
   else if ( north_ < 0.0 && south_ < 0.0) north_diff_south = fabs(north_) - fabs(south_);
   else north_diff_south  = fabs(north_) + fabs(south_);

   return ( north_diff_south/(rows() + 1) + 0.5);
}

long GribRegularLatLonGrid::computeIncLon() const
{
   return ((east_ - west_)/cols() + 0.5 );
}

long GribRegularLatLonGrid::computeRows(double north, double south, double west, double east) const
{
   if (north > 0.0 && south > 0.0 ) return (north - south)/nsIncrement_ + 1;
   else if ( north < 0.0 && south < 0.0) return (fabs(north) - fabs(south))/nsIncrement_ + 1;

   return (fabs(north) + fabs(south))/nsIncrement_ + 1;
}

long GribRegularLatLonGrid::computeCols(double west, double east) const
{
   return fabs((east - west)/weIncrement_) + 1;
}

// ========================================================================================
static PointList* read_number_of_data_points(grib_handle *h)
{
   // points to read
   long nb_nodes = 0;
   grib_get_long(h,"numberOfDataPoints",&nb_nodes);

   /// It should be noted that grib iterator is *only* available for certain grids
   /// i.e for Spherical Harmonics it is not implemented.
   int err = 0;
   grib_iterator *i = grib_iterator_new(h, 0, &err);
   if(err != 0 )
      throw std::runtime_error("error reading grib, could not create grib_iterator_new");

   PointList* pts = new PointList(nb_nodes);

   double lat   = 0.;
   double lon   = 0.;
   double value = 0.;

   size_t idx = 0;
   while( grib_iterator_next(i,&lat,&lon,&value) )
   {
      (*pts)[idx].assign(lat,lon);
      ++idx;
   }
   grib_iterator_delete(i);

   ASSERT( idx == nb_nodes );
   return pts;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
