/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_PolarStereoGraphic_H
#define atlas_grid_PolarStereoGraphic_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

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

   /// To compute bounding box in spherical co-ordinates we need:
   ///   1/ project the first spherical on to the plane.
   ///   2/ From this first x/y we compute the last point on the projected plane
   ///   3/ Convert the last point back into spherical space (lat/long)
   virtual BoundBox boundingBox() const;
   virtual size_t nPoints() const;

   /// The Point are in spherical (lat/long) co-ordinates
   /// Hence we need
   ///   1/ project the first spherical on to the plane. x/y
   ///   2/ Add x_grid_length_/y_grid_length_ to this point
   ///   3/ Convert this back to lat long values, and then repeat.
   virtual void coordinates( std::vector<double>& ) const;
   virtual void coordinates( std::vector<Point>& ) const;

   virtual std::string gridType() const;
   virtual GridSpec spec() const;
   virtual bool same(const Grid&) const;

private: // methods

   Point first_grid_pt() const { return first_grid_pt_;}
   long pts_along_xaxis() const { return npts_xaxis_;}
   long pts_along_yaxis() const { return npts_yaxis_;}
   long x_grid_length() const { return x_grid_length_; }
   long y_grid_length() const { return y_grid_length_; }


private: // members

   std::string hash_;

   Point first_grid_pt_;                  // This is in spherical lat long co-ordinate
   long npts_xaxis_;                     // No of points in x-axes *ON* the projected plane
   long npts_yaxis_;                     // No of points in y-axes *ON* the projected plane
   long x_grid_length_;                  // x grid length *ON* the projected plane
   long y_grid_length_;                  // y grid length *ON* the projected plane
   long lad_;                            // latitude where points npts_xaxis_ and npts_yaxis_ are specified
   long lov_;                            // orientation of grid
   bool north_pole_on_projection_plane_; // false means on south_pole_projection_plane
   bool spherical_earth_ ;               // true 6367470m,false mean oblate spheroid, 6378160m,6356775m, f=1/297.0
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
