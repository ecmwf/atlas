/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_PeriodicTransform_h
#define atlas_util_PeriodicTransform_h

#include "atlas/util/LonLatMicroDeg.h"

namespace atlas {
namespace util {

class PeriodicTransform
{
protected:
  double x_translation_;

public:
  PeriodicTransform()
  {
    x_translation_ = 360.;
  }

  void operator()(double source[2], double dest[2], double direction, double scale = 1.) const
  {
    dest[0] = source[0] + direction*x_translation_*scale;
    dest[1] = source[1];
  }

  void operator()(int source[2], int dest[2], int direction, int scale = 1) const
  {
    dest[0] = source[0] + direction*static_cast<int>(x_translation_*scale);
    dest[1] = source[1];
  }

  void operator()(double inplace[2], double direction, double scale = 1.) const
  {
    inplace[0] = inplace[0] + direction*x_translation_*scale;
    // inplace[1] = inplace[1]; null operation
  }

  void operator()(int inplace[2], int direction, int scale = 1) const
  {
    inplace[0] = inplace[0] + direction*static_cast<int>(x_translation_*scale);
    // inplace[1] = inplace[1]; null operation
  }

  void operator()(LonLatMicroDeg& inplace, int direction) const
  {
    inplace.set_lon( inplace.lon() + direction*microdeg(x_translation_) );
    // inplace.set_lat( inplace.lat() ); null operation
  }

};


} // namespace util
} // namespace atlas

#endif
