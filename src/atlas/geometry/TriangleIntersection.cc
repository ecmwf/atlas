/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "atlas/geometry/TriangleIntersection.h"

#include "eckit/eckit_config.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"

//----------------------------------------------------------------------------------------------------------------------

using Eigen::Vector3d;

namespace atlas {
namespace geometry {

Intersect TriangleIntersection::intersects(const Ray& r, double epsilon) const {

  Intersect isect;

  Vector3d edge1 = v1 - v0;
  Vector3d edge2 = v2 - v0;
  Vector3d pvec = r.dir.cross(edge2);

  const double det = edge1.dot(pvec);

  // ray is parallel to triangle (check ?)
  if (fabs(det) < epsilon) return isect.success(false);

  const double invDet = 1. / det;
  Vector3d tvec = r.orig - v0;
  isect.u = tvec.dot(pvec) * invDet;

  if(fabs(isect.u) < parametricEpsilon ) isect.u = 0;
  if(fabs(1-isect.u) < parametricEpsilon ) isect.u = 1;

  if (isect.u < 0 || isect.u > 1) return isect.success(false);

  Vector3d qvec = tvec.cross(edge1);
  isect.v = r.dir.dot(qvec) * invDet;

  if(fabs(isect.v) < parametricEpsilon ) isect.v = 0;
  if(fabs(1-isect.v) < parametricEpsilon ) isect.v = 1;

  double tmp = isect.u + isect.v;
  if(fabs(tmp) < parametricEpsilon ) tmp = 0;
  if(fabs(1-tmp) < parametricEpsilon ) tmp = 1;
  isect.u = tmp - isect.v;

  if (isect.v < 0 || isect.u + isect.v > 1) return isect.success(false);
  isect.t = edge2.dot(qvec) * invDet;

  return isect.success(true);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN
