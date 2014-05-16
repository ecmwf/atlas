#ifndef atlas_regular_lat_lon_grid_H
#define atlas_regular_lat_lon_grid_H
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
#include "grib_api.h"

#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class RegularLatLonGrid : public Grid {
public:
   RegularLatLonGrid( grib_handle* h );
   virtual ~RegularLatLonGrid();

   /// Overridden functions
   virtual std::string hash() const;
   virtual const char* gridType() const { return "regular_ll"; }
   virtual BoundBox boundingBox() const;
   virtual size_t nPoints() const { return points_.size(); }
   virtual void coordinates( Grid::Coords & ) const;
   /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
   virtual const std::vector<Point>& coordinates() const { return points_; }

   /// Functions specific to Regular Lat long grids
   long rows() const { return nptsNS_;}
   long cols() const { return nptsWE_;}
   Point latLon(size_t lat, size_t lon) const;
   double incLat() const { return nsIncrement_; }
   double incLon() const { return weIncrement_; }

private:
   // for verification/checks
   long computeIncLat() const ;
   long computeIncLon() const ;
   long computeRows(double north, double south, double west, double east) const;
   long computeCols(double west, double east) const;
   double epsilon() const;

private:
   double nsIncrement_;             /// In degrees
   double weIncrement_;             /// In degrees
   double north_;                   /// In degrees
   double south_;                   /// In degrees
   double west_;                    /// In degrees
   double east_;                    /// In degrees
   long nptsNS_;
   long nptsWE_;
   long   editionNumber_;           /// Grib 1 or Grib 2
   std::vector< Point > points_;     ///< storage of coordinate points
   std::string hash_;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
