/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grids_ReducedLonLatGrid_h
#define atlas_grids_ReducedLonLatGrid_h

#include "atlas/grids/ReducedGrid.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

/// @brief Reduced LonLat Grid
///
/// This grid is a special case of the class ReducedGrid, with
/// equidistant distribution of latitudes, and a equidistant distribution in zonal
/// direction, which reduce in number going closer towards poles,
/// essentially making the grid more uniform on the sphere
/// It can be constructed with following definition:
///   N   = number of latitudes in hemisphere
///   npts_per_lat[] = number of points on each latitude

class ReducedLonLatGrid: public ReducedGrid {

public:
  enum {EXCLUDES_POLES=0, INCLUDES_POLES=1};

private:

  struct defaults {
    // By default LonLat grids have the pole excluded
    static bool poles() { return EXCLUDES_POLES; }
  };

public:

  static std::string grid_type_str();

  ReducedLonLatGrid();

  ReducedLonLatGrid( const eckit::Parametrisation& );

  ReducedLonLatGrid( const size_t nlat, const long npts_per_lat[], bool poles = defaults::poles(), const Domain& domain = Domain::makeGlobal() );

  static std::string className();

  virtual eckit::Properties spec() const;

protected:

  void setup( const eckit::Parametrisation& );
  void setup( const size_t N, const long npts_per_lat[], bool poles=defaults::poles() );
  void set_typeinfo();

private:

  bool poles_;

};

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif // atlas_grids_ReducedLonLatGrid_h
