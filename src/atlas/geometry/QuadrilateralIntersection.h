/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_QuadIntersection_h
#define atlas_QuadIntersection_h

#include <limits>

#include "eckit/eckit_config.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"

#include "eckit/types/FloatCompare.h"

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

class QuadrilateralIntersection {

public:

  QuadrilateralIntersection(double* x0, double* x1, double* x2, double* x3) {

    v00 = Eigen::Vector3d::Map(x0);
    v10 = Eigen::Vector3d::Map(x1);
    v11 = Eigen::Vector3d::Map(x2);
    v01 = Eigen::Vector3d::Map(x3);

  }

  /// Intersection data structure

  struct Intersection {

    double u;
    double v;
    double t;

    void print(std::ostream& s) const { s << "Intersect[" << "(u=" << u << ",v=" << v << ",t=" << t << "]"; }

    friend std::ostream& operator<<(std::ostream& s, const Intersection& p) {
      p.print(s);
      return s;
    }

  };

  /// Ray trace data structure

  struct Ray {
    Eigen::Vector3d orig;
    Eigen::Vector3d dir;

    /// initializes ray with origin in point and direction to (0,0,0)
    explicit Ray(double* p) {
      orig = Eigen::Vector3d::Map(p);
      dir = -orig;
    }

    Ray(double* o, double* d) {
      orig = Eigen::Vector3d::Map(o);
      dir = Eigen::Vector3d::Map(d);
    }

    Eigen::Vector3d operator()(double t) const { return orig + t * dir; }
  };

  bool intersects(const Ray& r, Intersection& isect, const double epsilon = 2 * std::numeric_limits<double>::epsilon()) const;

private: // members

  Eigen::Vector3d v00;
  Eigen::Vector3d v10;
  Eigen::Vector3d v11;
  Eigen::Vector3d v01;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN

#endif
