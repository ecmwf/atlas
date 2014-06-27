/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_rotated_lat_lon_grid_H
#define atlas_grid_rotated_lat_lon_grid_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridFactory.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class RotatedLatLonGrid : public Grid {
   REGISTER(RotatedLatLonGrid);
public:
   RotatedLatLonGrid();
   virtual ~RotatedLatLonGrid();

   /// Overridden functions
   virtual std::string hash() const { return hash_;}
   virtual BoundBox boundingBox() const { return bbox_;}
   virtual size_t nPoints() const { return points_.size(); }
   virtual void coordinates( Grid::Coords & ) const;
   virtual std::string gridType() const { return std::string("rotated_ll"); }
   virtual GridSpec* spec() const;
   virtual void constructFrom(const GridSpec& );
   virtual bool compare(const Grid&) const;

   /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
   virtual const std::vector<Point>& coordinates() const { return points_; }

   /// Functions specific to Rotated Lat long grids
   double rotated_latitude() const { return rotated_latitude_; }
   double rotated_longitude() const { return rotated_longitude_; }
   double rotated_angle() const { return rotated_angle_; }

   Point latLon(size_t lat, size_t lon) const;
   long rows() const { return nptsNS_;}
   long cols() const { return nptsWE_;}
   double incLat() const { return nsIncrement_; }
   double incLon() const { return weIncrement_; }

private:
   std::string hash_;
   BoundBox bbox_;
   double rotated_latitude_;
   double rotated_longitude_;
   double rotated_angle_;

   double nsIncrement_;             ///< In degrees
   double weIncrement_;             ///< In degrees
   long nptsNS_;
   long nptsWE_;

   std::vector< Point > points_;     ///< storage of coordinate points

   /// Added friend mechanism to minimise data copying, during construction
   friend class GribRotatedLatLonGrid;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
