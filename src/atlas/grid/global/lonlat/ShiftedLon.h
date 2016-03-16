/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grids_global_lonlat_ShiftedLon_h
#define atlas_grids_global_lonlat_ShiftedLon_h

#include "atlas/grid/global/lonlat/ReducedLonLatGrid.h"

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------

class ShiftedLon: public ReducedLonLatGrid {

public:

  static std::string grid_type_str();

  ShiftedLon();

  /// @brief Constructor
  ShiftedLon( const eckit::Parametrisation& );

  /// @brief Constructor
  ShiftedLon( const size_t N );

  static std::string className();

  virtual eckit::Properties spec() const;

  size_t nlon() const { return ReducedGrid::nlon(0); }

  double lon( const size_t jlon ) const { return ReducedGrid::lon(0,jlon); }

protected:

  void setup( const eckit::Parametrisation& p);
  void setup( const size_t N );
  void set_typeinfo();
};

//------------------------------------------------------------------------------

} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas

#endif // atlas_grids_global_lonlat_ShiftedLon_h
