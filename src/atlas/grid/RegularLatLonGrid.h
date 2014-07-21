/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#ifndef atlas_grid_regular_lat_lon_grid_H
#define atlas_grid_regular_lat_lon_grid_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridFactory.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class RegularLatLonGrid : public Grid {
   REGISTER(RegularLatLonGrid);

public: // methods

   RegularLatLonGrid();
   virtual ~RegularLatLonGrid();

   /// Overridden functions
   virtual std::string hash() const { return hash_;}
   virtual BoundBox boundingBox() const { return bbox_;}
   virtual size_t nPoints() const { return points_.size(); }
   virtual void coordinates( Grid::Coords & ) const;
   virtual std::string gridType() const { return std::string("regular_ll"); }
   virtual GridSpec* spec() const;
   virtual void constructFrom(const GridSpec& );
   virtual void constructFrom(const eckit::Params&);
   virtual bool compare(const Grid&) const;

   /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
   virtual const std::vector<Point>& coordinates() const { return points_; }

   /// Functions specific to Regular Lat long grids
   Point latLon(size_t lat, size_t lon) const;
   long rows() const { return nptsNS_;}
   long cols() const { return nptsWE_;}

   double incLat() const { return nsIncrement_; }
   double incLon() const { return weIncrement_; }

   double computeIncLat() const;

   double computeIncLon() const;

   long computeRows(double north, double south, double west, double east) const;

   long computeCols(double west, double east) const;

   void computeCoords( std::vector< Point >& points_ );

private: // members

   std::string hash_;
   BoundBox bbox_;
   double nsIncrement_;             ///< In degrees
   double weIncrement_;             ///< In degrees
   long nptsNS_;
   long nptsWE_;
   std::vector< Point > points_;     ///< storage of coordinate points

   /// Added friend mechanism to minimise data copying, during construction
   friend class GribRegularLatLonGrid;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
