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

#include <vector>

#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

// gribs use the following convention: (from Shahram)
//
// Horizontally:  Points scan in the +i (+x) direction
// Vertically:    Points scan in the -j (-y) direction
//
// The way I verified this was to look at our SAMPLE files (which IFS uses).
// I also verified that IFS does not modify the scanning modes
// so whatever the samples say, is the convention
//
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

class RegularLatLon : public Grid {

public: // methods

	static std::string className() { return "atlas.grid.RegularLatLon"; }
	static std::string gridTypeStr() { return "regular_ll"; }

	RegularLatLon( const eckit::Params& p );

	RegularLatLon( size_t ni, size_t nj, const BoundBox& bbox );

	virtual ~RegularLatLon();

	virtual std::string uid() const;
	virtual std::string hash() const { return hash_;}

	virtual BoundBox boundingBox() const;
	virtual size_t nPoints() const;
	virtual void coordinates( std::vector<double>& ) const;
	virtual void coordinates( std::vector<Point>& ) const;

	virtual std::string gridType() const;
	virtual GridSpec* spec() const;

	virtual bool same(const Grid&) const;

protected: // methods

	long rows() const { return nptsNS_;}
	long cols() const { return nptsWE_;}

	double incLat() const { return incNS_; }
	double incLon() const { return incWE_; }

	double computeIncLat() const;

	double computeIncLon() const;

//	long computeRows(double north, double south, double west, double east) const;
//	long computeCols(double west, double east) const;

private: // members

	std::string hash_;

	BoundBox bbox_;

	double incNS_;             ///< in degrees
	double incWE_;             ///< in degrees

	long nptsNS_;
	long nptsWE_;

};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
