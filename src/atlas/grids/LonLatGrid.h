/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef atlas_grids_LonLatGrid_h
#define atlas_grids_LonLatGrid_h

#include "atlas/grids/ReducedLonLatGrid.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

/// @brief Regular LonLat Grid
///
/// This grid is a special case of the class ReducedGrid, with
/// equidistant distribution of latitudes and longitudes.
/// Longitude is the X-direction (first  index in C)
/// Latitude  is the Y-direction (second index in C)
class LonLatGrid: public ReducedLonLatGrid {

public:
  enum TYPE {EXCLUDES_POLES=0, INCLUDES_POLES=1};

private:

  struct defaults {
    // By default LonLat grids have the pole excluded
    static TYPE poles() { return EXCLUDES_POLES; }
  };

public:

  static std::string grid_type_str();

  LonLatGrid();

  /// @brief Constructor, possibly limited area
  ///
  /// If bounding box present in Params:
  ///   If (nlon,nlat) present:
  ///     dlon = (east-west)/(nlon-1)
  ///     dlat = (north-south)/(nlat-1)
  ///   Else
  ///     dlon = "lon_inc"
  ///     dlat = "lat_inc"
  ///   Longitudes: [west  :  dlon : east ]
  ///   Latitudes:  [north : -dlat : south]
  /// Else: global grid
  ///   - LonLatGrid(nlon,nlat,poles)
  ///   - LonLatGrid(londeg,latdeg,poles)
  LonLatGrid( const eckit::Parametrisation& );

  /// @brief Constructor, limited area grid
  ///
  /// dlon = (east-west)/(nlon-1)
  /// dlat = (north-south)/(nlat-1)
  /// Longitudes: [west  :  dlon : east ]
  /// Latitudes:  [north : -dlat : south]
  LonLatGrid( const size_t nlon, const size_t nlat, const BoundBox& );

  /// @brief Constructor, global grid
  ///
  /// If poles==false:
  ///   dlon = 360/nlon
  ///   dlat = 180/nlat
  ///   Longitudes: [0         :  dlon :  360-dlon ]
  ///   Latitudes:  [90-dlat/2 : -dlat : -90+dlat/2]
  /// Else:
  ///   dlon = 360/nlon
  ///   dlat = 180/(nlat-1)
  ///   Longitudes: [0  :  dlon :  360-dlon ]
  ///   Latitudes:  [90 : -dlat : -90       ]
  LonLatGrid( const size_t nlon, const size_t nlat, TYPE poles=defaults::poles() );

  /// @brief Constructor, global grid
  ///
  /// If poles==false:
  ///   Longitudes: [0           :  londeg :  360-londeg ]
  ///   Latitudes:  [90-latdeg/2 : -latdeg : -90+latdeg/2]
  /// Else:
  ///   Longitudes: [0  :  londeg :  360-londeg ]
  ///   Latitudes:  [90 : -latdeg : -90         ]
  LonLatGrid( const double londeg, const double latdeg, TYPE poles=defaults::poles() );

  /// @brief Constructor, limited area grid
  LonLatGrid( const double londeg, const double latdeg, const BoundBox& );

  /// @brief Constructor, global grid
  ///
  /// nlon = 2*nlat
  LonLatGrid( const size_t nlat, TYPE poles=defaults::poles() );

  static std::string className();

  virtual eckit::Properties spec() const;

  size_t nlon() const { return ReducedGrid::nlon(0); }

  double lon( const size_t jlon ) const { return ReducedGrid::lon(0,jlon); }

protected:

  void setup( const eckit::Parametrisation& p);
  void setup( const double londeg, const double latdeg, bool poles );
  void setup( const double londeg, const double latdeg, const BoundBox& );
  void setup( const size_t nlon, const size_t nlat, const BoundBox& );
  void setup( const size_t nlon, const size_t nlat, bool poles );
  void set_typeinfo();
};

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif // atlas_grids_LonLatGrid_h
