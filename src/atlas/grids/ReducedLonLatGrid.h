/*
 * (C) Copyright 1996-2014 ECMWF.
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
#include "atlas/Util.h"

namespace atlas {
namespace grids {

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

  static std::string gtype();

  ReducedLonLatGrid();

  ReducedLonLatGrid( const eckit::Params& );

  ReducedLonLatGrid( const int nlat, const int npts_per_lat[], bool poles=defaults::poles() );

  static std::string className();

  virtual GridSpec spec() const;

protected:

  void setup( const eckit::Params& );
  void setup( const int N, const int npts_per_lat[], bool poles=defaults::poles() );
  void set_typeinfo();

private:
  bool poles_;
};

register_BuilderT1(Grid,ReducedLonLatGrid,ReducedLonLatGrid::gtype());

} // namespace grids
} // namespace atlas

#endif // atlas_grids_ReducedLonLatGrid_h
