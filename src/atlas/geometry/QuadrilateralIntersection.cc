/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/geometry/QuadrilateralIntersection.h"

#include <cmath>

#include "eckit/eckit_config.h"
#include "eckit/exception/Exceptions.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"
#include "eckit/log/Log.h"

#include "atlas/geometry/TriangleIntersection.h"

using eckit::Log;
using Eigen::Vector3d;

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

double sign(const double& x)
{
  if( x >= 0. )
    return 1.0;
  else
    return -1.0;
}

/// @note Algorithm:
///       A. Lagae, P. Dutre', An Efficient Ray-Quadrilateral Intersection Test
///       JGTOOLS, 2005 vol. 10 (4) pp. 23-32.,
///       http://dx.doi.org/10.1080/2151237X.2005.10129208

/// @note variable names are maintained as in the pseudo-code block in the article

Intersect QuadrilateralIntersection::intersects(const Ray& r, double epsilon) const {

  Intersect isect; // intersection is false

  const Vector3d& O = r.orig;
  const Vector3d& D = r.dir;

  // rejects rays using the barycentric coords of intersection point wrt T

  Vector3d E01 = v10 - v00;
  Vector3d E03 = v01 - v00;

  Vector3d P = D.cross(E03);

  double det = E01.dot(P);

  if(fabs(det) < epsilon) return isect.success(false);

  Vector3d T = O - v00;

  double alpha = T.dot(P) / det;

  if(alpha < 0.) return isect.success(false);
  if(alpha > 1.) return isect.success(false);

  Vector3d Q = T.cross(E01);

  double beta = D.dot(Q) / det;

  if(beta < 0.) return isect.success(false);
  if(beta > 1.) return isect.success(false);

  // rejects rays using the barycentric coords of intersection point wrt T'
  // we resue the ' (prime) variables

  if((alpha+beta) > 1.) {

    Vector3d E23 = v01 - v11;
    Vector3d E21 = v10 - v11;

    Vector3d P = r.dir.cross(E21);

    double det = E23.dot(P);

    if(fabs(det) < epsilon) return isect.success(false);

    Vector3d T = O - v11;

    double alpha = T.dot(P) / det;

    if(alpha < 0.) return isect.success(false);

    Vector3d Q = T.cross(E23);

    double beta = D.dot(Q) / det;

    if(beta < 0.) return isect.success(false);
  }

  // compute the ray parameter of the intersection point

  isect.t = E03.dot(Q) / det;

  if(isect.t < 0) return isect.success(false);

  // compute the barycentric coordinates of V11

  Vector3d E02 = v11 - v00;

  Vector3d N = E01.cross(E03);

#define X 0
#define Y 1
#define Z 2

  const double Nx = fabs(N[X]);
  const double Ny = fabs(N[Y]);
  const double Nz = fabs(N[Z]);

  double alpha11 = 0.;
  double beta11  = 0.;

  if(Nx >= Ny && Nx >= Nz) {
    alpha11 = (E02[Y]*E03[Z] - E02[Z]*E03[Y]) / Nx;
    beta11  = (E01[Y]*E02[Z] - E01[Z]*E02[Y]) / Nx;
  }
  else
    if(Ny >= Nx && Ny >= Nz) {
      alpha11 = (E02[Z]*E03[X] - E02[X]*E03[Z]) / Ny;
      beta11  = (E01[Z]*E02[X] - E01[X]*E02[Z]) / Ny;
    }
    else {
      alpha11 = (E02[X]*E03[Y] - E02[Y]*E03[X]) / Nz;
      beta11  = (E01[X]*E02[Y] - E01[Y]*E02[X]) / Nz;
    }

#undef X
#undef Y
#undef Z

  // compute the bilinear coordinates of the intersection point

  if(fabs(alpha11 - 1.) < epsilon) {
    isect.u = alpha;
    if(fabs(beta11 - 1.) < epsilon) isect.v = beta;
    else isect.v = beta/(isect.u*(beta11-1.)+1.0);
  }
  else
    if(fabs(beta11 - 1.) < epsilon) {
        isect.v = beta;
        isect.u = alpha/(isect.v*(alpha11-1.)+1.);
    }
    else {
      double A = -(beta11 - 1.);
      double B = alpha*(beta11-1.) - beta*(alpha11-1.) - 1.0;
      double C = alpha;
      double Dt = B*B - 4*A*C;
      double Q = -0.5 * (B + sign(B)*sqrt(Dt));
      isect.u = Q/A;
      if(isect.u < 0. || isect.u > 1.) isect.u = C/Q;
      isect.v = beta/(isect.u*(beta11 - 1.) + 1.0);
    }

  return isect.success(true);
}

Intersect QuadrilateralIntersection::intersectsTG(const Ray &r, double epsilon) const
{
    Intersect isect; // intersection is false

    TriangleIntersection T013(v00.data(), v10.data(), v01.data());
    isect = T013.intersects(r,epsilon);
    if(isect)
        return isect;

    TriangleIntersection T231(v11.data(), v01.data(), v10.data());
    isect = T231.intersects(r,epsilon);
    if(isect)
    {
        isect.u = 1 - isect.u;
        isect.v = 1 - isect.v;
        return isect;
    }

    return isect.success(false);
}

bool QuadrilateralIntersection::validateIntersection(const Ray& r) const {

    TriangleIntersection T231(v11.data(), v01.data(), v10.data());
    if( T231.intersects(r) ) return true;

    TriangleIntersection T013(v00.data(), v10.data(), v01.data());
    if( T013.intersects(r) ) return true;

    TriangleIntersection T120(v10.data(), v11.data(), v00.data());
    if( T120.intersects(r) ) return true;

    TriangleIntersection T302(v01.data(), v00.data(), v11.data());
    if( T302.intersects(r) ) return true;

    return false;
}

bool QuadrilateralIntersection::validate() const {

    // normal for sub-triangle T231

    Vector3d E23 = v01 - v11;
    Vector3d E21 = v10 - v11;

    Vector3d N231 = E23.cross(E21);

    // normal for sub-triangle T013

    Vector3d E01 = v10 - v00;
    Vector3d E03 = v01 - v00;

    Vector3d N013 = E01.cross(E03);

    // normal for sub-triangle T120

    Vector3d E12 = - E21;
    Vector3d E10 = - E01;

    Vector3d N120 = E12.cross(E10);

    // normal for sub-triangle T302

    Vector3d E30 = - E03;
    Vector3d E32 = - E23;

    Vector3d N302 = E30.cross(E32);

    // all normals must point same way

    double dot02 = N231.dot(N013);
    double dot23 = N013.dot(N120);
    double dot31 = N120.dot(N302);
    double dot10 = N302.dot(N231);

//    Log::info() << dot02 << " "
//                << dot23 << " "
//                << dot31 << " "
//                << dot10 << std::endl;

    // all normals must point same way

    if (( dot02 >= 0. &&  dot23 >= 0. &&  dot31 >= 0. &&  dot10 >= 0. ) ||
        ( dot02 <= 0. &&  dot23 <= 0. &&  dot31 <= 0. &&  dot10 <= 0. ))
        return true;
    else
        return false;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN
