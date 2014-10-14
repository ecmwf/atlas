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
	static std::string gridTypeStr() { return "reduced_ll"; }

	ReducedLatLon( const eckit::Params& p );

	virtual ~ReducedLatLon();

	virtual std::string uid() const;
	virtual std::string hash() const { return hash_;}

	virtual BoundBox boundingBox() const { return bbox_;}
	virtual size_t nPoints() const;

	virtual void coordinates( std::vector<double>& ) const;
	virtual void coordinates( std::vector<Point>& ) const;

	virtual std::string gridType() const;
	virtual GridSpec spec() const;
	virtual bool same(const Grid&) const;

   void computeNPtsPerLat( std::vector<long>& );

protected: // methods

	long rows() const { return nptsNS_;}
	double incLat() const { return nsIncrement_; }

	long computeRows() const;


private: // members

	std::string hash_;
	BoundBox bbox_;
	double nsIncrement_;                   ///< In degrees
	long nptsNS_;                          ///< No of points along Y axes
	std::vector<long> nbPtsPerLat_;        ///< No of points per latitude
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
