/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_geometry_Triag3D_h
#define atlas_geometry_Triag3D_h

#include <limits>

#include "eckit/eckit_config.h"

#ifdef ECKIT_HAVE_EIGEN

#include "eckit/maths/Eigen.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/geometry/Ray.h"
#include "atlas/geometry/Intersect.h"

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

/// Triangle structure
/// Implements @link
/// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm

class Triag3D {

public: // types

  Triag3D(const double* x0, const double* x1, const double* x2) {

    v0 = Eigen::Vector3d::Map(x0);
    v1 = Eigen::Vector3d::Map(x1);
    v2 = Eigen::Vector3d::Map(x2);

  }

  Intersect intersects(const Ray& r, double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

  double area() const;

  void print(std::ostream& s) const { s << "Triag3D["
                                        << "v0="  << v0
                                        << ",v1=" << v1
                                        << ",v2=" << v2
                                        << "]"; }

  friend std::ostream& operator<<(std::ostream& s, const Triag3D& p) {
    p.print(s);
    return s;
  }

private: // members

  Eigen::Vector3d v0;
  Eigen::Vector3d v1;
  Eigen::Vector3d v2;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // ECKIT_HAVE_EIGEN

#endif
