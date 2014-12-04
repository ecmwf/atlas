/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_RotatedLatLon_H
#define atlas_RotatedLatLon_H

#include <cstddef>
#include <vector>

#include "eckit/memory/Builder.h"
#include "atlas/Grid.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

/// RotatedLatLon is a grid where the poles are shifted
///
///=== WMO specification ===
/// (6) Three parameters define a general latitude/longitude coordinate system,
///    formed by a general rotation of the sphere. One
///     choice for these parameters is:
///     (a)  The geographic latitude in degrees of the southern pole of the coordinate system, θp for example;
///
///     (b)  The geographic longitude in degrees of the southern pole of the coordinate system, λp for example;
///
///     (c)  The angle of rotation in degrees about the new polar axis
///          (measured clockwise when looking from the southern to the northern pole)
///          of the coordinate system, assuming the new
///          axis to have been obtained by first rotating the
///          sphere through λp degrees about the geographic polar axis, and
///          then rotating through (90 + θp) degrees so that
///          the southern pole moved along the (previously rotated) Greenwich meridian.
///=== end WMO specification ===

/// gribs use the following convention: (from Shahram)
///
/// Horizontally:  Points scan in the +i (+x) direction
/// Vertically:    Points scan in the -j (-y) direction
///
/// The way I verified this was to look at our SAMPLE files (which IFS uses).
/// I also verified that IFS does not modify the scanning modes
/// so whatever the samples say, is the convention
///
/// @todo Do we check the area? Can we assume area is multiple of the grids ?

class RotatedLatLon : public Grid {

public: // methods

	static std::string className() { return "atlas.grid.RotatedLatLon"; }
	static std::string grid_type_str() { return "rotated_ll"; }

	RotatedLatLon( const eckit::Params& p );
	virtual ~RotatedLatLon();

	virtual std::string uid() const;
	virtual std::string hash() const { return hash_;}

	virtual BoundBox bounding_box() const { return bbox_;}
	virtual size_t npts() const;

	virtual void lonlat( double[] ) const;
  virtual void lonlat( std::vector<double>& v ) const { Grid::lonlat(v); }
	virtual void lonlat( std::vector<Point>& ) const;

	virtual std::string grid_type() const;
	virtual GridSpec spec() const;
	virtual bool same(const Grid&) const;

private: // methods

	double rotated_latitude() const { return south_pole_lat_; }
	double rotated_longitude() const { return south_pole_lon_; }
	double rotated_angle() const { return south_pole_rot_angle_; }

	Point lonlat(size_t jlon, size_t jlat) const;
	long rows() const { return nptsNS_;}
	long cols() const { return nptsWE_;}
	double incLat() const { return nsIncrement_; }
	double incLon() const { return weIncrement_; }

private: // members

	std::string hash_;
	BoundBox bbox_;
	double south_pole_lat_;
	double south_pole_lon_;
	double south_pole_rot_angle_;

	double nsIncrement_;             ///< In degrees
	double weIncrement_;             ///< In degrees
	long nptsNS_;
	long nptsWE_;
};

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif
