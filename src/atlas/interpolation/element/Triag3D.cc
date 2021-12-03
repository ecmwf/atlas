/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>

#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/method/Intersect.h"
#include "atlas/interpolation/method/Ray.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace interpolation {
namespace element {

method::Intersect Triag3D::intersects(const method::Ray& r, double edgeEpsilon, double epsilon) const {
    method::Intersect isect;

    Vector3D edge1 = v1 - v0;
    Vector3D edge2 = v2 - v0;
    Vector3D pvec  = r.dir.cross(edge2);

    // ray is parallel to triangle (check?)
    const double det = edge1.dot(pvec);
    if (fabs(det) < epsilon) {
        return isect.fail();
    }

    const double invDet = 1. / det;

    Vector3D tvec = r.orig - v0;
    Vector3D qvec = tvec.cross(edge1);

    isect.u = tvec.dot(pvec) * invDet;
    isect.v = r.dir.dot(qvec) * invDet;
    isect.t = edge2.dot(qvec) * invDet;

    const double w = 1 - (isect.u + isect.v);

    if (w < 0) {
        // check if far outside of triangle, in respect to diagonal edge
        if (w < -edgeEpsilon) {
            return isect.fail();
        }

        // snap to diagonal edge
        // Note: we may still be to the left of vertical edge or below the
        // horizontal edge
        isect.u += 0.5 * w;
        isect.v += 0.5 * w;
    }

    if (isect.u < 0) {
        // check if far outside of triangle, in respect to vertical edge
        if ((isect.u < -edgeEpsilon) || (isect.v < -edgeEpsilon) || (isect.v > 1 + edgeEpsilon)) {
            return isect.fail();
        }

        // snap to lower/upper left corners
        isect.u = 0;
        if (isect.v < 0) {
            isect.v = 0;
        }
        else if (isect.v > 1) {
            isect.v = 1;
        }
    }

    if (isect.v < 0) {
        // check if far outside of triangle, in respect to horizontal edge
        if ((isect.v < -edgeEpsilon) || (isect.u < -edgeEpsilon) || (isect.u > 1 + edgeEpsilon)) {
            return isect.fail();
        }

        // snap to lower left/right corners
        isect.v = 0;
        if (isect.u < 0) {
            isect.u = 0;
        }
        else if (isect.u > 1) {
            isect.u = 1;
        }
    }

    return isect.success();
}

double Triag3D::area() const {
    Vector3D edge1 = v1 - v0;
    Vector3D edge2 = v2 - v0;

    Vector3D cross = edge1.cross(edge2);

    return 0.5 * cross.norm();
}

void Triag3D::print(std::ostream& s) const {
    s << "Triag3D["
      << "v0=(" << v0[0] << ", " << v0[1] << ", " << v0[2] << "), v1=(" << v1[0] << ", " << v1[1] << ", " << v1[2]
      << "), v2=(" << v2[0] << ", " << v2[1] << ", " << v2[2] << ")]";
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
