/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#ifndef atlas_regular_gaussian_grid_H
#define atlas_regular_gaussian_grid_H

#include <cstddef>
#include <vector>

#include "atlas/Grid.h"


//-----------------------------------------------------------------------------

namespace atlas {


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

	virtual std::string uid() const;
	virtual std::string hash() const { return hash_;}

	virtual BoundBox bounding_box() const { return bbox_; }
	virtual size_t npts() const { return npts_; }

	virtual void lonlat( double[] ) const;
	virtual void lonlat( std::vector<Point>& ) const;

	virtual std::string grid_type() const;
	virtual GridSpec spec() const;
	virtual bool same(const Grid&) const;

protected: // methods

	long gaussianNumber() const { return gaussN_;}

	long nj() const { return 2*gaussN_; }

	double computeIncLon() const;

	void computePoints( const std::vector<double>&, std::vector<Point>& pts ) const;
	long computeNPoints( const std::vector<double>& ) const;

	void computeLatitudes(std::vector<double>&) const;

private: // members

	std::string hash_;

	long        npts_;             ///< no of data points in grid, taking into account the bounding box
	long        gaussN_ ;                  ///< nb of points between pole and equator
	long        ni_;                       ///< nb of points along a parallel

	BoundBox    bbox_;

};

//-----------------------------------------------------------------------------


} // namespace eckit

#endif
