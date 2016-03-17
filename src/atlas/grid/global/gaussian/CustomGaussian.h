/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grids_global_gaussian_CustomGaussian_h
#define atlas_grids_global_gaussian_CustomGaussian_h

#include "atlas/grid/global/gaussian/Gaussian.h"

namespace atlas {
namespace grid {
namespace global {
namespace gaussian {

//------------------------------------------------------------------------------------------------------

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

class CustomGaussian: public Gaussian {
public:

  static std::string grid_type_str() { return "custom_gaussian"; }

  CustomGaussian( const eckit::Parametrisation& );

  CustomGaussian( const size_t N, const long npts_per_lat[] );

  static std::string className();

  virtual eckit::Properties spec() const;

protected:

  CustomGaussian() : Gaussian() {}

  void setup( const eckit::Parametrisation& );
  void setup_N_hemisphere( const size_t N, const long npts_per_lat[] );
  virtual void set_typeinfo();

};

//------------------------------------------------------------------------------------------------------

} // namespace gaussian
} // namespace global
} // namespace grid
} // namespace atlas

#endif // atlas_grids_global_gaussian_CustomGaussian_h
