/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#ifndef atlas_reduced_lat_lon_grid_H
#define atlas_reduced_lat_lon_grid_H

#include <cstddef>
#include <vector>

#include "atlas/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grids {

//-----------------------------------------------------------------------------

// ==================================================================================
// In a reduced gaussian grid, the number of longitude points along a latitude is specified.
// Latitudes may have differing numbers of points but the grid is symmetrical about the Equator.
// A reduced gaussian grid may also be called a quasi-regular gaussian grid.
//
// In the reduced grids used by ECMWF, the number of points on each latitude row
// is chosen so that the local east-west grid length remains approximately constant
// for all latitudes, with the restriction that the number should be suitable for the
// Fast Fourier Transform used to interpolate spectral fields to grid point fields,
// ie number = 2^p * 3^q * 5^r.
//
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

class ReducedLatLon : public Grid {

public: // methods

	static std::string className()   { return "atlas.grid.ReducedLatLon"; }
	static std::string gridTypeStr() { return "reduced_ll_depr"; }

	ReducedLatLon( const eckit::Params& p );

	virtual ~ReducedLatLon();

	virtual std::string uid() const;
	virtual std::string hash() const { return hash_;}

	virtual BoundBox bounding_box() const { return bbox_;}
	virtual size_t npts() const;

	virtual void lonlat( double[] ) const;
	virtual void lonlat( std::vector<Point>& ) const;

	virtual std::string grid_type() const;
	virtual GridSpec spec() const;
	virtual bool same(const Grid&) const;

private:
   void computenpts_per_lat( std::vector<int>& );

private: // methods

	long rows() const { return nptsNS_;}
	double incLat() const { return nsInc_; }

	long computeRows() const;
	size_t computeNPts() const;

	template <class T>
	void iterate(T&) const;

private: // members

	size_t                   npts_;        ///< no of data points in grid, taking into account the bounding box
	long                     nptsNS_;      ///< No of points along NS axis

	double                   nsInc_;       ///< increment in NS axis

	BoundBox                 bbox_;        ///< bounding box for data, only points within are considered part of grid

	std::vector<int> nbPtsPerLat_;        ///< No of points per latitude

	std::string hash_;
};

register_BuilderT1(Grid,ReducedLatLon,ReducedLatLon::gridTypeStr());

//-----------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif
