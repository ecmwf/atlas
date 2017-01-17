/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>

#include "atlas/interpolation/Ray.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace interpolation {

Ray::Ray(const double *p) {
    orig = Vector3D::Map(p);
    dir = -orig;
}

Ray::Ray(const double* o, const double* d) {
    orig = Vector3D::Map(o);
    dir = Vector3D::Map(d);
}

void Ray::print(std::ostream& s) const { s << "Ray[orig=" << orig << ",dir=" << dir << "]"; }

std::ostream& operator<<(std::ostream& s, const Ray& p) {
    p.print(s);
    return s;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace atlas

