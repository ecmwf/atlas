/*
 * (C) Copyright 1996-2016 ECMWF.
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

#include "atlas/geometry/Vector3D.h"
#include "atlas/geometry/Intersect.h"

namespace atlas {
namespace geometry {

class Ray;

//----------------------------------------------------------------------------------------------------------------------

class Quad3D {
public:

  Quad3D(const double* x0, const double* x1, const double* x2, const double* x3) {
    v00 = Vector3D::Map(x0);
    v10 = Vector3D::Map(x1);
    v11 = Vector3D::Map(x2);
    v01 = Vector3D::Map(x3);
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
  Vector3D v00; // aka v0
  Vector3D v10; // aka v1
  Vector3D v11; // aka v2
  Vector3D v01; // aka v3

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas


#endif
