/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#ifndef atlas_grid_reduced_lat_lon_grid_H
#define atlas_grid_reduced_lat_lon_grid_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class ReducedLatLonGrid : public Grid {
public:
   ReducedLatLonGrid();
   virtual ~ReducedLatLonGrid();

   /// Overridden functions
   virtual std::string hash() const { return hash_;}
   virtual BoundBox boundingBox() const { return bbox_;}
   virtual size_t nPoints() const { return points_.size(); }
   virtual void coordinates( Grid::Coords & ) const;
   virtual std::string gridType() const { return std::string("reduced_ll") ;}
   virtual GridSpec* spec() const;
   virtual void constructFrom(const GridSpec& );
   /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
   virtual const std::vector<Point>& coordinates() const { return points_; }

   /// Functions specific to Reduced Lat long grids
   long rows() const { return nptsNS_;}
   double incLat() const { return nsIncrement_; }

private:
   std::string hash_;
   BoundBox bbox_;
   double nsIncrement_;                   ///< In degrees
   long nptsNS_;                          ///< No of points along Y axes
   std::vector<long>    rgSpec_;          ///< No of points per latitude
   std::vector< Point > points_;           ///< storage of coordinate points

   /// Added friend mechanism to minimise data copying, during construction
   friend class GribReducedLatLonGrid;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
