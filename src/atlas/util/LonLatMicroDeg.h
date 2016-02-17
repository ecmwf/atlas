/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_LonLatMicroDeg_h
#define atlas_util_LonLatMicroDeg_h

#include "atlas/atlas_config.h"
#include "atlas/Parameters.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/Functions.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace util {

// ------------------------------------------------------------------------------------

// Forward declarations
class LonLatMicroDeg;
uidx_t unique_lonlat( const LonLatMicroDeg& );

// ------------------------------------------------------------------------------------

/// @brief LonLatMicroDegrees
///
/// Data is stored and accessed in microdegrees in int type (32 bits per coordinate component)
/// This structure is used to compare 2 lon-lat points in sorting algorithms
/// or to use in sending in MPI buffers using minimal memory,
/// The maximum range in degrees is [-2147.48 2147.48]
class LonLatMicroDeg
{
public:

// -- Constructors taking already microdegrees
  LonLatMicroDeg( int lonlat[2] )    { p[LON]=lonlat[LON]; p[LAT]=lonlat[LAT]; }
  LonLatMicroDeg( int lon, int lat ) { p[LON]=lon;         p[LAT]=lat; }

// -- Constructors taking degrees
  LonLatMicroDeg( const double& lon, const double& lat )  { p[LON]=microdeg(lon);         p[LAT]=microdeg(lat); }
  LonLatMicroDeg( const double lonlat[2] )                { p[LON]=microdeg(lonlat[LON]); p[LAT]=microdeg(lonlat[LAT]); }
  LonLatMicroDeg( const ArrayView<double,1>& lonlat )     { p[LON]=microdeg(lonlat[LON]); p[LAT]=microdeg(lonlat[LAT]); }
  LonLatMicroDeg( const eckit::geometry::Point2& lonlat ) { p[LON]=microdeg(lonlat[LON]); p[LAT]=microdeg(lonlat[LAT]); }

// -- Methods
  int lon() const { return p[LON]; }
  int lat() const { return p[LAT]; }
  void set_lon(int lon) { p[LON]=lon; }
  void set_lat(int lat) { p[LAT]=lat; }
  const int* data() const { return p; }
  uidx_t unique() { return unique_lonlat(*this); }

// -- Comparison operator
  bool operator < (const LonLatMicroDeg& other) const;

private:
  int p[2];
};

// ------------------------------------------------------------------------------------

inline bool LonLatMicroDeg::operator < (const LonLatMicroDeg& other) const
{
  if( p[LAT] > other.p[LAT]  ) return true;
  if( p[LAT] == other.p[LAT] ) return (p[LON] < other.p[LON]);
  return false;
}

// ------------------------------------------------------------------------------------

} // namespace util
} // namespace atlas

#include "atlas/util/Unique.h"

#endif
