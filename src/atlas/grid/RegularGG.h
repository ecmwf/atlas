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
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

class RegularGG : public Grid {

public: // methods

	static std::string className() { return "atlas.grid.RegularGG"; }
	static std::string gridTypeStr() { return "regular_gg"; }

	RegularGG( const eckit::Params& p );
	virtual ~RegularGG();

	virtual std::string hash() const { return hash_;}
	virtual BoundBox boundingBox() const { return bbox_; }
	virtual size_t nPoints() const { return points_.size(); }

	virtual void coordinates( std::vector<double>& ) const;
	virtual void coordinates( std::vector<Point>& ) const;

	virtual std::string gridType() const;
	virtual GridSpec* spec() const;
	virtual bool same(const Grid&) const;

protected: // methods

	Grid::Point latLon(size_t the_i, size_t the_j) const;
	long gaussianNumber() const { return gaussianNumber_;}

private: // members

	std::string hash_;
	BoundBox bbox_;
	std::vector< Point > points_;     ///< storage of coordinate points
	std::vector<double> latitudes_;
	long   gaussianNumber_;          /// No of points between pole and equator
	long   nj_;                     ///< No of points along Y axes

};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
