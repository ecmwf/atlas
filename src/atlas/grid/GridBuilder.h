#ifndef atlas_grid_builder_H
#define atlas_grid_builder_H
/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstddef>
#include <vector>
#include <cmath>

#include "grib_api.h"

#include "eckit/memory/NonCopyable.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/LocalPathName.h"

#include "eckit/geometry/Point2.h"
#include "atlas/grid/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

// Abstract base class for building Grid* objects
class GridBuilder : private eckit::NonCopyable {
public: // methods

   GridBuilder();
   virtual ~GridBuilder();

   virtual Grid::Ptr build(eckit::LocalPathName pathname) const = 0;
   // virtual Grid* build(eckit::DataHandle) = 0;
};

// =============================================================================

class GribGridBuilder : public GridBuilder {
private:
   GribGridBuilder();
public: // methods

   virtual ~GribGridBuilder();

   virtual Grid::Ptr build(eckit::LocalPathName pathname) const;


   Grid::Ptr build_grid_from_grib_handle(grib_handle* h ) const;

   /// @returns the singleton instance of this class
   static GribGridBuilder& instance();

   /// Returns the list of all known grib grid types.
   /// Allow better error handling during Grid construction.
   /// This really belongs in grib_api,
   static void known_grid_types(std::set<std::string>&);
};

// =============================================================================

class GribGridMaker : private eckit::NonCopyable {
public:
   virtual ~GribGridMaker();

   std::string hash() const { return hash_;}
   virtual Grid::BoundBox boundingBox() const = 0;
   virtual void coordinates( Grid::Coords& ) const = 0;
   virtual const std::vector<Grid::Point>& coordinates() const = 0;

protected:
   GribGridMaker( grib_handle* h );

   // grib edition 1 - milli-degrees
   // grib edition 2 - micro-degrees or could be defined by the keys: "subdivisionsOfBasicAngle" and "basicAngleOfTheInitialProductionDomain"
   // Therefore the client needs access to this when dealing with double based comparisons (for tolerances)
   double epsilon() const { return epsilon_; }

   static int scanningMode(long iScansNegatively, long jScansPositively);
   static void comparePointList(const std::vector<Grid::Point>& points,  double epsilon, grib_handle* handle);

protected:
   long   editionNumber_;           /// Grib 1 or Grib 2
   double north_;                   /// In degrees
   double south_;                   /// In degrees
   double west_;                    /// In degrees
   double east_;                    /// In degrees
   double epsilon_;                 /// Grib 1 or Grib 2
   long   numberOfDataPoints_;      ///
   long   iScansNegatively_;
   long   jScansPositively_;
   std::string hash_;
};

// =============================================================================

class GribReducedGaussianGrid : public GribGridMaker {
public:
   GribReducedGaussianGrid( grib_handle* h );
   virtual ~GribReducedGaussianGrid();

   /// Overridden functions
   virtual Grid::BoundBox boundingBox() const;
   virtual void coordinates( Grid::Coords & ) const;
   virtual const std::vector<Grid::Point>& coordinates() const { return points_; }

   /// Functions specific to Reduced Guassian grids
   long gaussianNumber() const { return gaussianNumber_;}
   const std::vector<double>& latitudes() { return latitudes_;}

private:
   void add_point(int lat_index);

   // TODO common code
   bool isGlobalNorthSouth() const;
   bool isGlobalWestEast() const;

private:
   long   gaussianNumber_;          /// No of points between pole and equator
   long   nj_;                      /// No of points along Y axes
   std::vector<long> rgSpec_;
   std::vector< Grid::Point > points_;     ///< storage of coordinate points
   std::vector<double> latitudes_;
};

// =============================================================================

class GribRegularGaussianGrid : public GribGridMaker {
public:
   GribRegularGaussianGrid( grib_handle* h );
   virtual ~GribRegularGaussianGrid();

   /// Overridden functions
   virtual Grid::BoundBox boundingBox() const;
   virtual void coordinates( Grid::Coords & ) const;
   virtual const std::vector<Grid::Point>& coordinates() const { return points_; }

   /// Functions specific to Regular Guassian grids
   Grid::Point latLon(size_t the_i, size_t the_j) const;
   long gaussianNumber() const { return gaussianNumber_;}
   const std::vector<double>& latitudes() { return latitudes_;}

private:
   // TODO these are common, common base class ?
   bool isGlobalNorthSouth() const;
   bool isGlobalWestEast() const;

private:
   long   nj_;                      /// No of points along Y axes
   long   gaussianNumber_;          /// No of points between pole and equator

   std::vector< Grid::Point > points_;     ///< storage of coordinate points
   std::vector<double> latitudes_;
};

// =============================================================================

class GribRegularLatLonGrid : public GribGridMaker {
public:
   GribRegularLatLonGrid( grib_handle* h );
   virtual ~GribRegularLatLonGrid();

   /// Overridden functions
   virtual Grid::BoundBox boundingBox() const;
   virtual size_t nPoints() const { return points_.size(); }
   virtual void coordinates( Grid::Coords & ) const;
   /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
   virtual const std::vector<Grid::Point>& coordinates() const { return points_; }

   /// Functions specific to Regular Lat long grids
   long rows() const { return nptsNS_;}
   long cols() const { return nptsWE_;}
   Grid::Point latLon(size_t lat, size_t lon) const;
   double incLat() const { return nsIncrement_; }
   double incLon() const { return weIncrement_; }

private:
   // for verification/checks
   long computeIncLat() const ;
   long computeIncLon() const ;
   long computeRows(double north, double south, double west, double east) const;
   long computeCols(double west, double east) const;

private:
   double nsIncrement_;             /// In degrees
   double weIncrement_;             /// In degrees
   long nptsNS_;
   long nptsWE_;
   std::vector< Grid::Point > points_;     ///< storage of coordinate points
};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
