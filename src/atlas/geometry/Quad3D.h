/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_geometry_Quad3D_h
#define atlas_geometry_Quad3D_h

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

class Quad3D {
public:

  Quad3D(const double* x0, const double* x1, const double* x2, const double* x3) {
    v00 = Eigen::Vector3d::Map(x0);
    v10 = Eigen::Vector3d::Map(x1);
    v11 = Eigen::Vector3d::Map(x2);
    v01 = Eigen::Vector3d::Map(x3);
  }

  Intersect intersects(const Ray& r, double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

  bool validate() const;

  double area() const;

  void print(std::ostream& s) const {
    s << "Quad3D[v00=" << v00 << ",v10=" << v10 << ",v11=" << v11 << ",v01=" << v01 << "]";
  }

  friend std::ostream& operator<<(std::ostream& s, const Quad3D& p) {
    p.print(s);
    return s;
  }

private:  // members
  Eigen::Vector3d v00; // aka v0
  Eigen::Vector3d v10; // aka v1
  Eigen::Vector3d v11; // aka v2
  Eigen::Vector3d v01; // aka v3

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // ECKIT_HAVE_EIGEN

#endif
