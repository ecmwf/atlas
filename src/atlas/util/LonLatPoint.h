/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_LonLatPoint_h
#define atlas_util_LonLatPoint_h

#include "atlas/atlas_config.h"
#include "atlas/Parameters.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/Functions.h"

namespace atlas {
namespace util {

struct LonLatPoint
{
  // Storage is in microdegrees
  // This structure is used in sorting algorithms, and uses less memory than
  // if x and y were in double precision.
  LonLatPoint() {}
  LonLatPoint( int x_, int y_ )
  {
    x = x_;
    y = y_;
  }
  LonLatPoint( long x_, long y_ )
  {
    x = x_;
    y = y_;
  }
  LonLatPoint( double x_, double y_ )
  {
    x = microdeg(x_);
    y = microdeg(y_);
  }
  LonLatPoint( const double coord[2] )
  {
    x = microdeg(coord[LON]);
    y = microdeg(coord[LAT]);
  }
  LonLatPoint( const ArrayView<int,1>& coord )
  {
    x = coord[LON];
    y = coord[LAT];
  }
  LonLatPoint( const ArrayView<double,1>& coord )
  {
    x = microdeg(coord[LON]);
    y = microdeg(coord[LAT]);
  }

  long uid64() const;

  int uid32() const;

  gidx_t uid() const;

  mutable int x, y;
  bool operator < (const LonLatPoint& other) const
  {
    if( y > other.y  ) return true;
    if( y == other.y ) return (x < other.x);
    return false;
  }
private:

  template<typename T> gidx_t uidT() const;

  static int WEST;  
  static int EAST;
  static int NORTH;
  static int SOUTH;
};

template<> inline gidx_t LonLatPoint::uidT<int >() const { return uid32(); }
template<> inline gidx_t LonLatPoint::uidT<long>() const { return uid64(); }

inline int LonLatPoint::uid32() const
{
  // max precision is 0.02 degree
  int iy = static_cast<int>((2*NORTH-y)*5e-5);
  int ix = static_cast<int>((x+2*EAST)*5e-5);
  iy <<= 17;
  int id = iy | ix;
  return id;
}

inline long LonLatPoint::uid64() const
{
  // max precision is 1 microdegree
  long iy = static_cast<long>((4*NORTH-y));
  long ix = static_cast<long>((x+4*EAST));
  iy <<= 31;
  long id = iy | ix;
  return id;
}

inline gidx_t LonLatPoint::uid() const {
  return uidT<gidx_t>();
}


} // namespace util
} // namespace atlas

#endif
