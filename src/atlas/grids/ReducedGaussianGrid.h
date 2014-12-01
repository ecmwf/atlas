/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_grids_ReducedGaussianGrid_h
#define atlas_grids_ReducedGaussianGrid_h

#include "atlas/grids/ReducedGrid.h"
#include "atlas/Util.h"

namespace atlas {
namespace grids {

/// @brief Reduced Gaussian Grid
///
/// This grid is a special case of the class ReducedGrid, in which
/// the latitudes are distributed according to the roots of the
/// Legendre Polynomials, and a equidistant distribution in zonal
/// direction, which reduce in number going closer towards poles,
/// essentially making the grid more uniform on the sphere
/// It can be constructed with following definition:
///   N   = number of latitudes in hemisphere
///   npts_per_lat[] = number of points on each latitude

class ReducedGaussianGrid: public ReducedGrid {

public:

  static std::string gtype() { return "reduced_gg"; }

  ReducedGaussianGrid();

  ReducedGaussianGrid( const eckit::Params& );

  ReducedGaussianGrid( const int N, const int npts_per_lat[] );

  static std::string className();

  virtual GridSpec spec() const;

protected:

  void setup( const eckit::Params& );
  void setup_N_hemisphere( const int N, const int npts_per_lat[] );
  void set_typeinfo();

};

register_BuilderT1(Grid,ReducedGaussianGrid,ReducedGaussianGrid::gtype());

} // namespace grids
} // namespace atlas

#endif // atlas_grids_ReducedGaussianGrid_h
