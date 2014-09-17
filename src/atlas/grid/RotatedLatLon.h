/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_RotatedLatLon_H
#define atlas_grid_RotatedLatLon_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"


//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
/// RotatedLatLon is a grid where the poles are shifted
///
/// @note
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
	static std::string gridTypeStr() { return "rotated_ll"; }

	RotatedLatLon( const eckit::Params& p );
	virtual ~RotatedLatLon();

	virtual std::string uid() const;
	virtual std::string hash() const { return hash_;}

	virtual BoundBox boundingBox() const { return bbox_;}
	virtual size_t nPoints() const;

	virtual void coordinates( std::vector<double>& ) const;
	virtual void coordinates( std::vector<Point>& ) const;

	virtual std::string gridType() const;
	virtual GridSpec spec() const;
	virtual bool same(const Grid&) const;

protected: // methods

	double rotated_latitude() const { return south_pole_lat_; }
	double rotated_longitude() const { return south_pole_lon_; }
	double rotated_angle() const { return south_pole_rot_angle_; }

	Point latLon(size_t lat, size_t lon) const;
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

//-----------------------------------------------------------------------------

//rotgrid
// See; http://proj.badc.rl.ac.uk/svn/cows/cows_support/rotated_grid_utils/trunk/lib/rotated_grid_utils/rotgrid.py
//
// Contains routines to calculate the position of a given point as
// seen in a rotated grid.
//
// The easiest interface is by instantiating Rotgrid and calling the
// transform() method.
//
// For repeated calculations with the same lats and lons (e.g. transforming
// every point on the non-rotated grid), it may be more efficient to call
// rotgrid_core directly.  This requires some terms to have been
// pre-calculated by the calling routine.

class Rotgrid {
public:
//    Rotated grid class.  For more info, see doc strings for '__init__'
//    and 'transform' methods.

   // Set up rotated grid for transformations.
   //        Inputs:
   //
   //        south_pole_lon, south_pole_lat: longitude (degrees) and latitude (degrees) of the
   //                          pole of the rotated grid, as seen in the
   //                          non-rotated grid
   //
   //        polerotate: optional input -- by default, the calculation assumes
   //                    that the rotated grid is singly rotated, i.e. that the
   //                    common meridian which passes through poles of rotated
   //                    and non-rotated grid has the same longitude value in
   //                    both grids.  If there is additional rotation about the
   //                    pole of the rotated grid, then set this input to the
   //                    value in degrees
   //
   //        nPoleGridLon: an alternative way of specifying the longitudinal
   //                      rotation between grids: specify as the longitude
   //                      (degrees) of the true north pole as seen in the
   //                      rotated grid.  If set, overrides polerotate.
   //
   //        lonMin:  minimum longitude for output of transforms to be perfomed
   //                 defaults to -180 so that longitudes are normally output in
   //                 the range [-180, 180) but e.g. specify as 0 if [0, 360) is
   //                 desired.
   Rotgrid(const Grid::Point& south_pole,
           double polerotate = 0,
           double nPoleGridLon = 0,
           double lonMin = 0); // for -180 - +180, choose lonMin as -180

   // Performs transformations to/from rotated grid.
   //
   //        Inputs:
   //
   //          lon, lat: longitude (degrees) and latitude (degrees)
   //                    of a point X, as seen in the non-rotated grid
   //
   //          inverse: optional input -- set to a true value for inverse transform
   //                          (coords on rotated grid to coords on nonrotated)
   //
   //        Returns:
   //
   //            The coordinates of the point X (in degrees) as seen in the rotated
   //            grid (or the non-rotated grid in case of inverse transform), as a
   //            2-element tuple: (longitude, latitude)
   Grid::Point transform(const Grid::Point& latlon, bool inverse = false) const;


   // magics, assume south_pole_rot_angle_ = 0;
   Grid::Point magics_rotate( const Grid::Point&) const;
   Grid::Point magics_unrotate( const Grid::Point&) const;

   Grid::Point rotate( const Grid::Point&) const;
   Grid::Point unrotate( const Grid::Point&) const;

private:
   // Inputs:
   //
   //      cossouth_pole_lat, sinsouth_pole_lat:
   //            cos and sine of latitude of the pole
   //            of the rotated grid, as seen in the non-rotated grid
   //
   //      sindlon, cosdlon, coslat, sinlat:
   //            cos and sine of longitude offset
   //            and cos and sine of latitude
   //            of a point X, as seen in the non-rotated grid
   //
   //            (NB longitude offset is taken from the common meridian which
   //            passes through poles of rotated and non-rotated grid)
   //
   //    Returns:
   //
   //      The coordinates of the point X (in radians) as seen in the rotated grid.
   //      as a 2-element tuple: (longitude offset, latitude)
   //
   //      (NB longitude offset is taken from the common meridian which
   //      passes through poles of the rotated of non-rotated grids)
   Grid::Point rotgrid_core(
               double cossouth_pole_lat, double sinsouth_pole_lat,
               double cosdlon, double sindlon, double coslat, double sinlat) const;

private:
   double degree_to_radian_;
   double radian_to_degree_;

   Grid::Point south_pole_;
   double south_pole_rot_angle_;

   double cossouth_pole_lat_;
   double sinsouth_pole_lat_;
   double lonmin_;
   double lonmax_;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
