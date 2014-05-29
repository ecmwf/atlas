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

#include "atlas/grid/GridBuilder.h"
#include "atlas/grid/Unstructured.h"

#include "eckit/log/Log.h"
#include "eckit/geometry/Point3.h"
#include "eckit/grib/GribAccessor.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/thread/Mutex.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/io/StdFile.h"
#include "eckit/exception/Exceptions.h"

using namespace std;
using namespace eckit;
using namespace eckit::geometry;


namespace atlas {
namespace grid {

typedef std::vector< atlas::grid::Grid::Point > PointList;
static PointList* read_number_of_data_points(grib_handle *h);

//=====================================================================================

GridBuilder::GridBuilder()  {}
GridBuilder::~GridBuilder() {}

//=====================================================================================

static Mutex local_mutex;
static Mutex mutex_instance;

GRIBGridBuilder::GRIBGridBuilder()  {}
GRIBGridBuilder::~GRIBGridBuilder() {}

Grid::Ptr GRIBGridBuilder::build(const eckit::PathName& pathname) const
{
   AutoLock<Mutex> lock(local_mutex);

   /// Will by default open in read mode, destructor will close file
   StdFile theFile(pathname);

   int err;
   grib_handle* handle = grib_handle_new_from_file(0,theFile,&err);
   if (err != 0 || handle == 0) {
      std::stringstream ss; ss << "Could not create grib handle for file " << pathname << "  err=" << err;
      throw SeriousBug(ss.str(),Here());
   }

   Grid::Ptr the_grid_ptr = build_grid_from_grib_handle(handle);

   err = grib_handle_delete(handle);
   if (err != 0) {
      std::stringstream ss; ss << "Could not delete handle for file " << pathname << " ";
      throw SeriousBug(ss.str(),Here());
   }

   return the_grid_ptr;
}


Grid::Ptr GRIBGridBuilder::build_grid_from_grib_handle( grib_handle* handle ) const
{
   char string_value[64];
   size_t len = sizeof(string_value)/sizeof(char);
   if (grib_get_string(handle,"gridType",string_value,&len)  != 0) {
      throw SeriousBug(string("grib_get_string failed for gridType"),Here());
   }

   if (strncasecmp(string_value,"regular_ll",10) == 0) {
      GribRegularLatLonGrid maker( handle );
      return maker.build();
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
      return maker.build();
   }
   else if (strncasecmp(string_value,"regular_gg",10) == 0) {
      GribRegularGaussianGrid maker( handle );
      return maker.build();
   }

   // Unknown grid type, extract data points from the grib handle
   PointList* pntlist = read_number_of_data_points(handle);
   return Grid::Ptr(new grid::Unstructured( pntlist, grib_hash(handle) ));
}


GRIBGridBuilder& GRIBGridBuilder::instance()
{
   AutoLock<Mutex> lock(mutex_instance);

   static GRIBGridBuilder* obj = 0;

   if( !obj )
      obj = new GRIBGridBuilder();

   return *obj;
}

void GRIBGridBuilder::known_grid_types(std::set<std::string>& grids)
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

GribGridBuilderHelper::GribGridBuilderHelper(grib_handle* handle)
: handle_(handle),
  editionNumber_(0),
  north_(0.0),south_(0.0),west_(0.0),east_(0.0),
  epsilon_(1e-6),
  numberOfDataPoints_(0)
{
   if (handle == NULL) {
      throw SeriousBug(string("NULL grib_handle"),Here());
   }

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

GribGridBuilderHelper::~GribGridBuilderHelper(){}


int GribGridBuilderHelper::scanningMode(long iScansNegatively, long jScansPositively)
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

void GribGridBuilderHelper::comparePointList(const std::vector<Grid::Point>& points,  double epsilon, grib_handle* handle)
{
   // Check point list compared with grib
   PointList* pointlist1 = read_number_of_data_points(handle);
   const std::vector<atlas::grid::Grid::Point>& grib_pntlist = *pointlist1;
   ASSERT( points.size() == grib_pntlist.size());

   int print_point_list = 0;
   std::streamsize old_precision = cout.precision();
   for(size_t i =0;  i < points.size() && i < grib_pntlist.size(); i++) {
      if (!FloatCompare::is_equal(points[i].lat(),grib_pntlist[i].lat(),epsilon) ||
          !FloatCompare::is_equal(points[i].lon(),grib_pntlist[i].lon(),epsilon))
      {
         if (print_point_list == 0) Log::info() << " Point list DIFFER, show first 10, epsilon(" << epsilon << ")" << std::endl;
         Log::info() << setw(3) << i << " :"
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << points[i].lat() << ", "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << points[i].lon() << "  "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << grib_pntlist[i].lat() << ", "
                  << setw(20) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << grib_pntlist[i].lon()
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
: GribGridBuilderHelper(handle),
  nj_(0),
  the_grid_( new ReducedGaussianGrid() )
{
   Log::info() << "Build a GribReducedGaussianGrid  " << std::endl;

   the_grid_->bbox_ = boundingBox();
   the_grid_->hash_ = hash_;
}

Grid::Ptr GribReducedGaussianGrid::build()
{
   // Extract the guassian grid attributes from the grib handle

   GRIB_CHECK(grib_get_long(handle_,"numberOfParallelsBetweenAPoleAndTheEquator",&the_grid_->gaussianNumber_),0);
   GRIB_CHECK(grib_get_long(handle_,"Nj",&nj_),0);

   // get reduced grid specification. These are number of point's along the lines of latitude
   size_t rgSpecLength = 0;
   GRIB_CHECK(grib_get_size(handle_,"pl",&rgSpecLength),0);
   the_grid_->rgSpec_.resize(rgSpecLength);
   GRIB_CHECK(grib_get_long_array(handle_,"pl",&(the_grid_->rgSpec_[0]),&rgSpecLength),0);
   ASSERT( rgSpecLength == nj_);


   // This provides the 'y' co-ordinate of each line of latitude
   the_grid_->latitudes_.resize(2 *the_grid_->gaussianNumber_ );
   grib_get_gaussian_latitudes(the_grid_->gaussianNumber_, &the_grid_->latitudes_[0]);

   // number of lines of latitude should, be twice the numberOfParallelsBetweenAPoleAndTheEquator
   ASSERT( the_grid_->rgSpec_.size() == 2 * the_grid_->gaussianNumber_);


   // Create point list based on area. To avoid rounding errors determine if we have
   // Global area, then use simple algorithm.
   if (isGlobalNorthSouth() && isGlobalWestEast()) {
      Log::info() << " GLOBAL " << std::endl;
      for ( int i = 0; i < 2*the_grid_->gaussianNumber_ && i < the_grid_->rgSpec_.size(); i++ ) {
         long no_of_points_along_latitude = the_grid_->rgSpec_[i];
         double east_west_grid_length = 360.0/no_of_points_along_latitude;
         double plon = 0;
         for(int k = 0; k < no_of_points_along_latitude; k++) {
            the_grid_->points_.push_back( Grid::Point( the_grid_->latitudes_[i], plon ) );
            plon += east_west_grid_length;
         }
      }
   }
   else {
      Log::info() << " LOCAL " << std::endl;
      for ( int i = 0; i < 2*the_grid_->gaussianNumber_ && i < the_grid_->rgSpec_.size(); i++ ) {
         if ( FloatCompare::is_equal(the_grid_->latitudes_[i],north_,epsilon()) ) {
            // Log::info() << " same north latitudes_[i] " <<  latitudes_[i] << " north_ " << north_ << std::endl;
            add_point(i);
            continue;
         }
         if ( FloatCompare::is_equal(the_grid_->latitudes_[i],south_,epsilon()) ) {
            add_point(i);
            continue;
         }
         if ( the_grid_->latitudes_[i] < north_ && the_grid_->latitudes_[i] > south_) {
            add_point(i);
         }
      }
   }
   if (jScansPositively_ == 1 )
      std::reverse(the_grid_->latitudes_.begin(), the_grid_->latitudes_.end());

   double EXPECTED_longitudeOfLastGridPointInDegrees = 360.0 - (90.0/(the_grid_->gaussianNumber_));
   Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
   Log::info() << " epsilon()                                      " << epsilon() << std::endl;
   Log::info() << " gaussianNumber_                                " << the_grid_->gaussianNumber_ << std::endl;
   Log::info() << " nj                                             " << nj_ << std::endl;
   Log::info() << " isGlobalNorthSouth()                           " << isGlobalNorthSouth() << std::endl;
   Log::info() << " isGlobalWestEast()                             " << isGlobalWestEast() << std::endl;
   Log::info() << " pl-size                                        " << rgSpecLength <<  " pl[0]="<< the_grid_->rgSpec_[0] << " pl[1]=" << the_grid_->rgSpec_[1] << std::endl;
   Log::info() << " latitudes_.size ()                             " << the_grid_->latitudes_.size() << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively_ << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively_ << std::endl;
   Log::info() << " scanning_mode                                  " << scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << east_ << std::endl;
   Log::info() << " EXPECTED longitudeOfLastGridPointInDegrees     " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << EXPECTED_longitudeOfLastGridPointInDegrees << std::endl;
   Log::info() << " numberOfDataPoints                             " << numberOfDataPoints_ << std::endl;
   int no_of_points_in_pl = 0;
   for(int i = 0; i < the_grid_->rgSpec_.size(); i++) no_of_points_in_pl += the_grid_->rgSpec_[i];
   Log::info() << " no_of_points_in_pl                             " << no_of_points_in_pl << std::endl;
   Log::info() << " points_.size()                                 " << the_grid_->points_.size() << std::endl;

   ASSERT(nj_ == 2*the_grid_->gaussianNumber_);
   ASSERT(the_grid_->points_.size() == numberOfDataPoints_);
   ASSERT(FloatCompare::is_equal(east_,EXPECTED_longitudeOfLastGridPointInDegrees,epsilon()));

   // Check point list compared with grib
   comparePointList(the_grid_->points_,epsilon(),handle_);

   return Grid::Ptr( the_grid_.release() );
}

GribReducedGaussianGrid::~GribReducedGaussianGrid()
{
   Log::info() << "Destroy a ReducedGaussianGrid" << std::endl;
}

void GribReducedGaussianGrid::add_point(int lat_index)
{
   long no_of_points_along_latitude = the_grid_->rgSpec_[lat_index];
   double east_west_grid_length = 360.0/no_of_points_along_latitude;
   double plon = 0;
   for(int k = 0; k < no_of_points_along_latitude; k++) {
      if ( (plon > west_ && plon < east_) ||
           (FloatCompare::is_equal(plon,west_,epsilon())  ||
            FloatCompare::is_equal(plon,east_,epsilon()) )
      ) {
         the_grid_->points_.push_back( Grid::Point( the_grid_->latitudes_[lat_index], plon ) );
         plon += east_west_grid_length;
      }
   }
}

Grid::BoundBox GribReducedGaussianGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox( Grid::Point(south_,west_),Grid::Point(north_,east_));
}

bool GribReducedGaussianGrid::isGlobalNorthSouth() const
{
   return (the_grid_->gaussianNumber_*2 == nj_);
}

bool GribReducedGaussianGrid::isGlobalWestEast() const
{
   // GRIB way of determining globalness
   if (west_ == 0) {
      double last_long = 360.0 - (90.0/(double)the_grid_->gaussianNumber_) ;
      return FloatCompare::is_equal(east_,last_long,epsilon());
   }
   return false;
}

// ========================================================================================

GribRegularGaussianGrid::GribRegularGaussianGrid(grib_handle* handle)
: GribGridBuilderHelper(handle),
  nj_(0),
  the_grid_( new RegularGaussianGrid() )
{
   Log::info() << "Build a RegularGaussianGrid  " << std::endl;

   the_grid_->bbox_ = boundingBox();
   the_grid_->hash_ = hash_;
}

GribRegularGaussianGrid::~GribRegularGaussianGrid()
{
   Log::info() << "Destroy a RegularGaussianGrid" << std::endl;
}

Grid::Ptr GribRegularGaussianGrid::build()
{
   // Extract the guassian grid attributes from the grib handle

    GRIB_CHECK(grib_get_long(handle_,"numberOfParallelsBetweenAPoleAndTheEquator",&the_grid_->gaussianNumber_),0);
    GRIB_CHECK(grib_get_long(handle_,"Nj",&nj_),0);

    double array[2 *the_grid_->gaussianNumber_];
    grib_get_gaussian_latitudes(the_grid_->gaussianNumber_, array);


    double nptsWE = 4 * the_grid_->gaussianNumber_ ;
    double weIncrement = 360.0/nptsWE;
    if (isGlobalNorthSouth() && isGlobalWestEast()) {
       Log::info() << " GLOBAL                      weIncrement = " << weIncrement << std::endl;
       for ( int i = 0; i < 2*the_grid_->gaussianNumber_; i++ ) {
          double lat = array[i];
          ASSERT(lat < 90.0 && lat > -90.0);
          the_grid_->latitudes_.push_back(lat);
       }
       if (jScansPositively_ == 1 )
          std::reverse(the_grid_->latitudes_.begin(), the_grid_->latitudes_.end());
       for(size_t i = 0 ; i < the_grid_->latitudes_.size(); i++) {
          double plon = 0;
          for( size_t j = 0; j < nptsWE; ++j) {
             the_grid_->points_.push_back( Grid::Point( the_grid_->latitudes_[i], plon ) );
             plon += weIncrement;
          }
       }
    }
    else {
       Log::info() << " LOCAL " << std::endl;
       for ( int i = 0; i < 2*the_grid_->gaussianNumber_; i++ ) {
          double lat = array[i];
          ASSERT(lat < 90.0 && lat > -90.0);
          if ( FloatCompare::is_equal(lat,north_,epsilon()) ) {
             the_grid_->latitudes_.push_back(lat);
             continue;
          }
          if ( FloatCompare::is_equal(lat,south_,epsilon()) ) {
             the_grid_->latitudes_.push_back(lat);
             continue;
          }
          if ( lat < north_ && lat > south_) {
             the_grid_->latitudes_.push_back(lat);
          }
       }

       if (jScansPositively_ == 1 )
          std::reverse(the_grid_->latitudes_.begin(), the_grid_->latitudes_.end());
       for(size_t i = 0 ; i < the_grid_->latitudes_.size(); i++) {
          double plon = west_;
          for( size_t j = 0; j < nptsWE; ++j) {
             if (plon < east_ || FloatCompare::is_equal(plon,east_)) {
                the_grid_->points_.push_back( Grid::Point( the_grid_->latitudes_[i], plon ) );
                plon += weIncrement;
             }
          }
       }
    }

    double EXPECTED_longitudeOfLastGridPointInDegrees = 360.0 - (90.0/(the_grid_->gaussianNumber_));
    Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
    Log::info() << " epsilon()                                      " << epsilon() << std::endl;
    Log::info() << " gaussianNumber_                                " << the_grid_->gaussianNumber_ << std::endl;
    Log::info() << " iScansNegatively                               " << iScansNegatively_ << std::endl;
    Log::info() << " jScansPositively                               " << jScansPositively_ << std::endl;
    Log::info() << " scanning_mode                                  " << scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
    Log::info() << " latitudeOfFirstGridPointInDegrees              " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << north_ << std::endl;
    Log::info() << " longitudeOfFirstGridPointInDegrees             " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << west_ << std::endl;
    Log::info() << " latitudeOfLastGridPointInDegrees               " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << south_ << std::endl;
    Log::info() << " longitudeOfLastGridPointInDegrees              " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << east_ << std::endl;
    Log::info() << " EXPECTED longitudeOfLastGridPointInDegrees     " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << EXPECTED_longitudeOfLastGridPointInDegrees << std::endl;
    Log::info() << " iDirectionIncrementInDegrees(west-east   incr) " << weIncrement << std::endl;
    Log::info() << " nptsNS                                         " << 2 * the_grid_->gaussianNumber_ << std::endl;
    Log::info() << " nptsWE                                         " << nptsWE << std::endl;
    Log::info() << " numberOfDataPoints                             " << numberOfDataPoints_ << std::endl;
    Log::info() << " -----------------------------------------------" << std::endl;
    Log::info() << " points_.size()  " << the_grid_->points_.size() << "       numberOfDataPoints " << numberOfDataPoints_ << std::endl;
    Log::info() << " point[0]                               " << the_grid_->points_[0].lat() << ", " << the_grid_->points_[0].lon() <<  std::endl;
    Log::info() << " point[1]                               " << the_grid_->points_[1].lat() << ", " << the_grid_->points_[1].lon() <<  std::endl;

    ASSERT(the_grid_->points_.size() == numberOfDataPoints_);

    ASSERT(FloatCompare::is_equal(east_,EXPECTED_longitudeOfLastGridPointInDegrees,epsilon()));

    // Check point list compared with grib
    comparePointList(the_grid_->points_,epsilon(),handle_);

    // take ownership
    return Grid::Ptr(the_grid_.release() );
}

Grid::BoundBox GribRegularGaussianGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Grid::Point(south_,west_),Grid::Point(north_,east_));
}

bool GribRegularGaussianGrid::isGlobalNorthSouth() const
{
   return (the_grid_->gaussianNumber_*2 == nj_);
}

bool GribRegularGaussianGrid::isGlobalWestEast() const
{
   // GRIB way of determining globalness
   if (west_ == 0) {
      double last_long = 360.0 - (90.0/(double)the_grid_->gaussianNumber_) ;
      return FloatCompare::is_equal(east_,last_long,epsilon());
   }
   return false;
}

// =====================================================================================

GribRegularLatLonGrid::GribRegularLatLonGrid(grib_handle* handle)
: GribGridBuilderHelper(handle),
  the_grid_( new RegularLatLonGrid() )
{
   Log::info() << "Build a RegularLatLonGrid  " << std::endl;
   the_grid_->bbox_ = boundingBox();
   the_grid_->hash_ = hash_;
}

GribRegularLatLonGrid::~GribRegularLatLonGrid()
{
   Log::info() << "Destroy a RegularLatLonGrid" << std::endl;
}

Grid::Ptr GribRegularLatLonGrid::build()
{
   // Extract the regular lat long grid attributes from the grib handle

   GRIB_CHECK(grib_get_double(handle_,"jDirectionIncrementInDegrees",&(the_grid_->nsIncrement_)),0);
   GRIB_CHECK(grib_get_double(handle_,"iDirectionIncrementInDegrees",&(the_grid_->weIncrement_)),0);

   GRIB_CHECK(grib_get_long(handle_,"Nj",&(the_grid_->nptsNS_)),0);
   GRIB_CHECK(grib_get_long(handle_,"Ni",&(the_grid_->nptsWE_)),0);


   double plat = north_;
   the_grid_->points_.reserve( (the_grid_->nptsNS_ + 1) * (the_grid_->nptsWE_ + 1) );
   for( size_t j = 0; j < the_grid_->nptsNS_; ++j) {
      double plon = west_;
      for( size_t i = 0; i < the_grid_->nptsWE_; ++i) {
         the_grid_->points_.push_back( Grid::Point( plat, plon ) );
         plon += the_grid_->weIncrement_;
      }
      plat -= the_grid_->nsIncrement_;
   }

   Log::info() << " editionNumber                                  " << editionNumber_ << std::endl;
   Log::info() << " iScansNegatively                               " << iScansNegatively_ << std::endl;
   Log::info() << " jScansPositively                               " << jScansPositively_ << std::endl;
   Log::info() << " scanning_mode                                  " << scanningMode(iScansNegatively_,jScansPositively_) << std::endl;
   Log::info() << " latitudeOfFirstGridPointInDegrees              " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << north_ << std::endl;
   Log::info() << " longitudeOfFirstGridPointInDegrees             " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << west_ << std::endl;
   Log::info() << " latitudeOfLastGridPointInDegrees               " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << south_ << std::endl;
   Log::info() << " longitudeOfLastGridPointInDegrees              " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << east_ << std::endl;
   Log::info() << " jDirectionIncrementInDegrees(north-south incr) " << the_grid_->nsIncrement_ << std::endl;
   Log::info() << " iDirectionIncrementInDegrees(west-east   incr) " << the_grid_->weIncrement_ << std::endl;
   Log::info() << " Nj(num of points North South)                  " << the_grid_->nptsNS_ << std::endl;
   Log::info() << " Ni(num of points West East)                    " << the_grid_->nptsWE_ << std::endl;
   Log::info() << " numberOfDataPoints                             " << numberOfDataPoints_ << std::endl;
   Log::info() << " -----------------------------------------------" << std::endl;
   Log::info() << " computeIncLat() " << computeIncLat() << "      nsIncrement_ " << the_grid_->nsIncrement_ << std::endl;
   Log::info() << " computeIncLon() " << computeIncLon() << "      weIncrement_ " << the_grid_->nsIncrement_ << std::endl;
   Log::info() << " computeRows()   " << computeRows(north_,south_,west_,east_) << "     nptsNS_ " << the_grid_->nptsNS_ << std::endl;
   Log::info() << " computeCols()   " << computeCols(west_,east_) <<  "     nptsWE_ " << the_grid_->nptsWE_ << std::endl;
   Log::info() << " points_.size()  " << the_grid_->points_.size() << "     numberOfDataPoints_ " << numberOfDataPoints_ << std::endl << std::endl;

   ASSERT(the_grid_->nsIncrement_ == computeIncLat());
   ASSERT(the_grid_->weIncrement_ == computeIncLon());
   ASSERT(the_grid_->nptsNS_ == computeRows(north_,south_,west_,east_));
   ASSERT(the_grid_->nptsWE_ == computeCols(west_,east_));
   ASSERT(the_grid_->points_.size() == numberOfDataPoints_);

   // Check point list compared with grib
   comparePointList(the_grid_->points_,epsilon(),handle_);

   // take ownership
   return Grid::Ptr( the_grid_.release() );
}

Grid::BoundBox GribRegularLatLonGrid::boundingBox() const
{
   // Grid::BoundBox expects bottom left, top right
   return Grid::BoundBox(Grid::Point(south_,west_),Grid::Point(north_,east_));
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
   if (north > 0.0 && south > 0.0 ) return (north - south)/the_grid_->nsIncrement_ + 1;
   else if ( north < 0.0 && south < 0.0) return (fabs(north) - fabs(south))/the_grid_->nsIncrement_ + 1;

   return (fabs(north) + fabs(south))/the_grid_->nsIncrement_ + 1;
}

long GribRegularLatLonGrid::computeCols(double west, double east) const
{
   return fabs((east - west)/the_grid_->weIncrement_) + 1;
}

long GribRegularLatLonGrid::rows() const { return the_grid_->rows();}
long GribRegularLatLonGrid::cols() const { return the_grid_->cols();}
double GribRegularLatLonGrid::incLat() const { return the_grid_->incLat(); }
double GribRegularLatLonGrid::incLon() const { return the_grid_->incLon(); }

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
   if(err != 0 ) {
      throw SeriousBug(string("Error reading grib. Could not create grib_iterator_new"),Here()) ;
   }

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
