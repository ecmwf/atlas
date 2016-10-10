/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <cmath>
#include "atlas/interpolation/Ray.h"
#include "atlas/interpolation/Triag3D.h"

//----------------------------------------------------------------------------------------------------------------------


static const double parametricEpsilon = 1e-7; ///< Epsilon used to compare weights and u,v's


namespace atlas {
namespace interpolation {

Intersect Triag3D::intersects(const Ray& r, double epsilon) const {

    Intersect isect;

    Vector3D edge1 = v1 - v0;
    Vector3D edge2 = v2 - v0;
    Vector3D pvec = r.dir.cross(edge2);

    // ray is parallel to triangle (check?)
    const double det = edge1.dot(pvec);
    if (fabs(det) < epsilon) return isect.success(false);

    // pick an epsilon based on a small proportion of the smallest edge length
    // (this scales linearly so it better compares with linear weights u,v,w)
    const double edgeEpsilon = std::max(epsilon, parametricEpsilon * std::sqrt(std::min(edge1.norm2(), edge2.norm2())));

    const double invDet = 1. / det;
    Vector3D tvec = r.orig - v0;
    isect.u = tvec.dot(pvec) * invDet;

    if (fabs(isect.u) < edgeEpsilon) isect.u = 0;
    else if (fabs(1-isect.u) < edgeEpsilon) isect.u = 1;
    else if (isect.u < 0 || isect.u > 1) return isect.success(false);

    Vector3D qvec = tvec.cross(edge1);
    isect.v = r.dir.dot(qvec) * invDet;

    if (fabs(isect.v) < edgeEpsilon) isect.v = 0;
    else if (fabs(1-isect.v) < edgeEpsilon) isect.v = 1;
    else if (isect.v < 0 || isect.v > 1) return isect.success(false);

    // check if projection point is near the diagonal of reference triangle
    // if very close (but outside) to diagonal, equally correct u and v to the inside
    const double w = 1 - (isect.u + isect.v);
    if (w < 0) {
        if (fabs(w) < edgeEpsilon) {
            isect.u += 0.5 * w;
            isect.v += 0.5 * w;
        } else {
            return isect.success(false);
        }
    }

    isect.t = edge2.dot(qvec) * invDet;

    return isect.success(true);
}

double Triag3D::area() const
{
    Vector3D edge1 = v1 - v0;
    Vector3D edge2 = v2 - v0;

    Vector3D cross = edge1.cross(edge2);

    return 0.5 * cross.norm();
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace atlas

