/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_PolarStereoGraphic_H
#define atlas_PolarStereoGraphic_H

#include <cstddef>
#include <vector>

#include "atlas/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grids {

//-----------------------------------------------------------------------------
/// PolarStereoGraphic is a projected grid
///

class PolarStereoGraphic : public Grid {

public: // methods

   static std::string className() { return "atlas.grid.PolarStereoGraphic"; }
   static std::string gridTypeStr() { return "polar_stereographic"; }

   PolarStereoGraphic( const eckit::Params& p );
   virtual ~PolarStereoGraphic();

   virtual std::string uid() const;
   virtual std::string hash() const { return hash_;}

   /// To compute bounding box in spherical co-ordinates we need to:
   ///   1/ project the first spherical pt, on to the plane.
   ///   2/ From this first x/y we compute the last point on the projected plane
   ///   3/ Convert the last point back into spherical space (lat/long)
   virtual BoundBox bounding_box() const;
   virtual size_t npts() const;

   /// The Points are in spherical (lat/long) co-ordinates
   /// Hence we need:
   ///   1/ project the first spherical point on to the plane x/y
   ///   2/ Add x_grid_length_/y_grid_length_ to this point
   ///   3/ Convert this back to lat long values, and then repeat.
   virtual void lonlat( double[] ) const;
   virtual void lonlat( std::vector<double>& v ) const { Grid::lonlat(v); }
   virtual void lonlat( std::vector<Point>& ) const;

   virtual std::string grid_type() const;
   virtual GridSpec spec() const;
   virtual bool same(const Grid&) const;

private: // members

   std::string hash_;

   bool iScansPositively_;            // Used to determine correct bounding box
   bool jScansPositively_;            // Used to determine correct bounding box
   Point first_grid_pt_;               // This is in spherical lat long co-ordinate
   long npts_xaxis_;                  // No of points in x-axes *ON* the projected plane
   long npts_yaxis_;                  // No of points in y-axes *ON* the projected plane
   long x_grid_length_;               // x grid length *ON* the projected plane, in meters
   long y_grid_length_;               // y grid length *ON* the projected plane, in meters
   long resolutionAndComponentFlag_;  // is used to determine sphere/oblate, But this is extracted separately,
                                       // to avoid having to mess with bits. Needed to match geography hash grid_spec -> grid in tests
   double lad_;                       // latitude where points npts_xaxis_ and npts_yaxis_ are specified
   double orientationOfTheGrid_;      // east longitude value, in degrees ( longitude of natural origin)
   bool southPoleOnProjectionPlane_;
   bool earth_is_oblate_;             // true 6367470m,false mean oblate spheroid, 6378160m,6356775m, f=1/297.0
   double radius_;                    // default 6371229
   double semi_major_;                // default 6378137   (a)
   double semi_minor_;                // default 6356752.3 (b)
   double e_;                         // calculated e = sqrt( 1 - b*b/a*a)
};

//-----------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif
