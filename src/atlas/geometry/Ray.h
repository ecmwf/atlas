/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_geometry_Ray_h
#define atlas_geometry_Ray_h

#include <limits>

#include "eckit/eckit_config.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"

#include "eckit/types/FloatCompare.h"

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

/// Ray trace data structure

struct Ray {

  Eigen::Vector3d orig;
  Eigen::Vector3d dir;

  /// initializes ray with origin in point and direction to (0,0,0)
  explicit Ray(const double* p);

  Ray(const double* o, const double* d) {
    orig = Eigen::Vector3d::Map(o);
    dir = Eigen::Vector3d::Map(d);
  }

  Eigen::Vector3d operator()(double t) const { return orig + t * dir; }

  void print(std::ostream& s) const { s << "Ray[orig=" << orig << ",dir=" << dir << "]"; }

  friend std::ostream& operator<<(std::ostream& s, const Ray& p) {
    p.print(s);
    return s;
  }

};

/// Intersection data structure

struct Intersect {

  double u;
  double v;
  double t;

  Intersect() : u(0.), v(0.), t(0.), success_(false) {}

  operator bool() const { return success_; }

  Intersect& success(bool s){ success_ = s; return *this; }

  void print(std::ostream& s) const { s << "Intersect[u=" << u << ",v=" << v << ",t=" << t <<",success=" << success_ << "]"; }

  friend std::ostream& operator<<(std::ostream& s, const Intersect& p) {
    p.print(s);
    return s;
  }

private:

  bool success_;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN

#endif
