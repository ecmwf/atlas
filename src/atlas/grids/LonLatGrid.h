/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_grids_LonLatGrid_h
#define atlas_grids_LonLatGrid_h

#include "atlas/grids/ReducedGrid.h"
#include "atlas/Util.h"

namespace atlas {
namespace grids {

/// @brief Regular LonLat Grid
///
/// This grid is a special case of the class ReducedGrid, with
/// equidistant distribution of latitudes and longitudes.
/// Longitude is the X-direction (first  index in C)
/// Latitude  is the Y-direction (second index in C)
class LonLatGrid: public ReducedGrid {

public:

  static std::string gtype();

  LonLatGrid();

  /// @brief Constructor, possibly limited area
  ///
  /// If bounding box present in Params:
  ///   dlon = (east-west)/(nlon-1)
  ///   dlat = (north-south)/(nlat-1)
  ///   Longitudes: [west  :  dlon : east ]
  ///   Latitudes:  [north : -dlat : south]
  /// Else:
  ///   As constructor LonLatGrid(nlon,nlat) for global grid
  ///
  /// The constructor with bounding box will *NOT* create a cropped version
  /// of the global one, but rather creates a new (nlon x nlat) grid within
  /// given bounding box... This is somewhat inconsistent with GaussianGrid
  /// behaviour, where the limited area fits exactly in the global grid.
  /// @todo introduce new parameter to make distinction.
  LonLatGrid( const eckit::Params& );

  /// @brief Constructor, limited area grid
  ///
  /// dlon = (east-west)/(nlon-1)
  /// dlat = (north-south)/(nlat-1)
  /// Longitudes: [west  :  dlon : east ]
  /// Latitudes:  [north : -dlat : south]
  LonLatGrid( const int nlon, const int nlat, const BoundBox& );

  /// @brief Constructor, global grid
  ///
  /// dlon = 360/nlon
  /// dlat = 180/nlat
  /// Longitudes: [0         :  dlon :  360-dlon ]
  /// Latitudes:  [90-dlat/2 : -dlat : -90+dlat/2]
  LonLatGrid( const int nlon, const int nlat );

  /// @brief Constructor, global grid
  ///
  /// N = nlat/2 = number of latitudes between pole and equator
  /// nlon = 4*N
  /// nlat = 2*N
  LonLatGrid( const int N );

  static std::string className();

  virtual GridSpec spec() const;

  int nlon() const { return ReducedGrid::nlon(0); }

  double lon( const int jlon ) const { return ReducedGrid::lon(0,jlon); }

protected:

  void setup(const eckit::Params& p);
  void setup(const int nlon, const int nlat, const BoundBox& );
  void setup(const int nlon, const int nlat);
  void set_typeinfo();
};

register_BuilderT1(Grid,LonLatGrid,LonLatGrid::gtype());

} // namespace grids
} // namespace atlas

#endif // atlas_grids_LonLatGrid_h
