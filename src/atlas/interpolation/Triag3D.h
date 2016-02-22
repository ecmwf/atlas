/*
 * (C) Copyright 1996-2016 ECMWF.
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

#include "atlas/interpolation/Vector3D.h"
#include "atlas/interpolation/Intersect.h"

namespace atlas {
namespace geometry {

class Ray;

//----------------------------------------------------------------------------------------------------------------------

/// Triangle structure
/// Implements @link
/// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm

class Triag3D {

public: // types

  Triag3D(const Vector3D& x0, const Vector3D& x1, const Vector3D& x2):
    v0(x0),
    v1(x1),
    v2(x2) {
  }

  Triag3D(const double* x0, const double* x1, const double* x2) {

    v0 = Vector3D::Map(x0);
    v1 = Vector3D::Map(x1);
    v2 = Vector3D::Map(x2);

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

  Vector3D v0;
  Vector3D v1;
  Vector3D v2;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas


#endif
