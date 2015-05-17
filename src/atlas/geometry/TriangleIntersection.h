/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_TriangleIntersection_h
#define atlas_TriangleIntersection_h

#include <limits>

#include "eckit/eckit_config.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"

#include "eckit/types/FloatCompare.h"


namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

/// Triangle structure

class TriangleIntersection {

public: // types

  TriangleIntersection(double* x0, double* x1, double* x2) {

    v0 = Eigen::Vector3d::Map(x0);
    v1 = Eigen::Vector3d::Map(x1);
    v2 = Eigen::Vector3d::Map(x2);

  }

  /// Intersection data structure

  struct Intersection {

    double u;
    double v;
    double t;

    double w() const { return 1.0 - u - v; }

    void print(std::ostream& s) const { s << "Intersect[u=" << u << ",v=" << v << ",w=" << w() << ",t=" << t << "]"; }

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

  /// Implements @link
  /// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm
  bool intersects(const Ray& r,
                  Intersection& isect,
                  const double epsilon = 2 * std::numeric_limits<double>::epsilon()) const;

private: // members

  Eigen::Vector3d v0;
  Eigen::Vector3d v1;
  Eigen::Vector3d v2;


};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN

#endif
