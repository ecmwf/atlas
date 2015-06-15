/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grids_GaussianGrid_h
#define atlas_grids_GaussianGrid_h

#include "atlas/grids/ReducedGaussianGrid.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

/// @brief Regular Gaussian Grid
///
/// This grid is a special case of the class ReducedGrid, in which
/// the latitudes are distributed according to the roots of the
/// Legendre Polynomials, and a equidistant distribution in zonal
/// direction, equal for every latitude.
/// It can be constructed with following definition:
///   N   = number of latitudes in hemisphere
///   4*N = number of longitudes along one latitude
///
/// A gaussian grid is a latitude/longitude grid.
/// The spacing of the latitudes is not regular.
/// However, the spacing of the lines of latitude is symmetrical about the Equator.
/// Note that there is no latitude at either Pole or at the Equator.
/// A grid is usually referred to by its 'number' N, which is the number of lines of latitude between a Pole and the Equator.
/// The longitudes of the grid points are defined by giving the number of points along each line of latitude.
/// The first point is at longitude 0 and the points are equally spaced along the line of latitude.
/// In a regular Gaussian grid, the number of longitude points along a latitude is 4*N.
/// In a reduced Gaussian grid, the number of longitude points along a latitude is specified.
/// Latitudes may have differing numbers of points but the grid is symmetrical about the Equator.
/// A reduced gaussian grid may also be called a quasi-regular Gaussian grid.
///
/// ==================================================================================
/// gribs use the following convention: (from Shahram)
///
/// Horizontally:  Points scan in the +i (+x) direction
/// Vertically:    Points scan in the -j (-y) direction
///
/// The way I verified this was to look at our SAMPLE files (which IFS uses).
/// I also verified that IFS does not modify the scanning modes
/// so whatever the samples say, is the convention
/// ==================================================================================
/// Area: Do we check the area.
/// Area: Can we assume area is multiple of the grids ?
class GaussianGrid: public ReducedGaussianGrid {

public:

  static std::string grid_type_str() { return "regular_gg"; }

  GaussianGrid();

  GaussianGrid( const eckit::Params& );

  GaussianGrid( const int N );

  static std::string className();

  virtual GridSpec spec() const;

  int nlon() const { return ReducedGaussianGrid::nlon(0); }

  double lon( const int jlon ) const { return ReducedGaussianGrid::lon(0,jlon); }

protected:

  void setup(const int N);
  void setup_lat_hemisphere(const int N, const double lats[]);
  void set_typeinfo();
};

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif // atlas_grids_GaussianGrid_h
