/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/geometry/TriangleIntersection.h"

#include "eckit/eckit_config.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"

//----------------------------------------------------------------------------------------------------------------------

using Eigen::Vector3d;

namespace atlas {
namespace geometry {

bool TriangleIntersection::intersects(const Ray& r, Intersection& isect, const double epsilon) const {

  Vector3d edge1 = v1 - v0;
  Vector3d edge2 = v2 - v0;
  Vector3d pvec = r.dir.cross(edge2);

  const double det = edge1.dot(pvec);

  if (std::abs(det) < epsilon) return false;

  const double invDet = 1. / det;
  Vector3d tvec = r.orig - v0;
  isect.u = tvec.dot(pvec) * invDet;
  if (isect.u + epsilon < 0 || isect.u - epsilon > 1) return false;

  Vector3d qvec = tvec.cross(edge1);
  isect.v = r.dir.dot(qvec) * invDet;

  if (isect.v + epsilon < 0 || isect.u + isect.v - epsilon > 1) return false;
  isect.t = edge2.dot(qvec) * invDet;

  return true;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN
