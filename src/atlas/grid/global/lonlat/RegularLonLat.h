/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grids_RegularLonLat_h
#define atlas_grids_RegularLonLat_h

#include "atlas/grid/global/lonlat/LonLat.h"

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------------------------------

/// @brief Regular LonLat Grid
///
/// This grid is a special case of the class ReducedGrid, with
/// equidistant distribution of latitudes and longitudes.
/// Longitude is the X-direction (first  index in C)
/// Latitude  is the Y-direction (second index in C)
class RegularLonLat: public ReducedLonLatGrid {

public:

  static std::string grid_type_str();

  /// @brief Constructor
  ///
  ///   - RegularLonLat(nlon,nlat)
  ///   - RegularLonLat(londeg,latdeg)
  RegularLonLat( const eckit::Parametrisation& );

  /// @brief Constructor, global grid
  ///
  ///   dlon = 360/nlon
  ///   dlat = 180/(nlat-1)
  ///   Longitudes: [0  :  dlon :  360-dlon ]
  ///   Latitudes:  [90 : -dlat : -90       ]
  RegularLonLat( const size_t nlon, const size_t nlat );

  /// @brief Constructor, global grid
  ///
  ///   Longitudes: [0  :  londeg :  360-londeg ]
  ///   Latitudes:  [90 : -latdeg : -90         ]
  RegularLonLat( const double &londeg, const double &latdeg );

  /// @brief Constructor, global grid
  ///
  /// nlon = 4*N
  /// nlat = 2*N+1
  /// londeg = latdeg = 90/N
  RegularLonLat( const size_t N );

  static std::string className();

  virtual eckit::Properties spec() const;

  size_t nlon() const { return ReducedGrid::nlon(0); }

  double lon( const size_t jlon ) const { return ReducedGrid::lon(0,jlon); }

protected:

  void setup( const eckit::Parametrisation& p);
  void setup( const double londeg, const double latdeg );
  void setup( const size_t nlon, const size_t nlat );
  void set_typeinfo();
};

//------------------------------------------------------------------------------------------------------

} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas

#endif // atlas_grids_RegularLonLat_h
