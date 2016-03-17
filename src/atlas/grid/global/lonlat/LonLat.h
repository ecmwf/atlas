/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grids_global_lonlat_LonLat_h
#define atlas_grids_global_lonlat_LonLat_h

#include "atlas/grid/global/Structured.h"

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------

/// @brief Reduced LonLat Grid
///
/// This grid is a special case of the class ReducedGrid, with
/// equidistant distribution of latitudes, and a equidistant distribution in zonal
/// direction, which reduce in number going closer towards poles,
/// essentially making the grid more uniform on the sphere
/// It can be constructed with following definition:
///   N   = number of latitudes in hemisphere
///   npts_per_lat[] = number of points on each latitude

class LonLat: public ReducedGrid {

public:

  static std::string grid_type_str();

  LonLat();

  static std::string className();

protected:

  virtual void set_typeinfo() = 0;

};

//------------------------------------------------------------------------------

} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas

#endif // atlas_grids_global_lonlat_LonLat_h
