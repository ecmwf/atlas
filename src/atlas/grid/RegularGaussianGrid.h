/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#ifndef atlas_grid_regular_gaussian_grid_H
#define atlas_grid_regular_gaussian_grid_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridFactory.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

// A gaussian grid is a latitude/longitude grid.
// The spacing of the latitudes is not regular.
// However, the spacing of the lines of latitude is symmetrical about the Equator.
// Note that there is no latitude at either Pole or at the Equator.
// A grid is usually referred to by its 'number' N, which is the number of lines of latitude between a Pole and the Equator.
// The longitudes of the grid points are defined by giving the number of points along each line of latitude.
// The first point is at longitude 0 and the points are equally spaced along the line of latitude.
// In a regular Gaussian grid, the number of longitude points along a latitude is 4*N.
// In a reduced Gaussian grid, the number of longitude points along a latitude is specified.
// Latitudes may have differing numbers of points but the grid is symmetrical about the Equator.
// A reduced gaussian grid may also be called a quasi-regular Gaussian grid.

class RegularGaussianGrid : public Grid {
   REGISTER(RegularGaussianGrid);
public:
   RegularGaussianGrid();
   virtual ~RegularGaussianGrid();

   /// Overridden functions
   virtual std::string hash() const { return hash_;}
   virtual BoundBox boundingBox() const { return bbox_; }
   virtual size_t nPoints() const { return points_.size(); }
   virtual void coordinates( Grid::Coords & ) const;
   virtual std::string gridType() const { return std::string("regular_gg"); }
   virtual GridSpec* spec() const;
   virtual void constructFrom(const GridSpec& );

   /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
   virtual const std::vector<Point>& coordinates() const { return points_; }

   /// Functions specific to Regular Guassian grids
   Grid::Point latLon(size_t the_i, size_t the_j) const;
   long gaussianNumber() const { return gaussianNumber_;}

private:
   std::string hash_;
   BoundBox bbox_;
   std::vector< Point > points_;     ///< storage of coordinate points
   std::vector<double> latitudes_;
   long   gaussianNumber_;          /// No of points between pole and equator
   long   nj_;                     ///< No of points along Y axes

   /// Added friend mechanism to minimise data copying, during construction
   friend class GribRegularGaussianGrid;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
