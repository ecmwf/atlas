/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <cmath>

#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/interpolation/method/Intersect.h"
#include "atlas/interpolation/method/Ray.h"

//----------------------------------------------------------------------------------------------------------------------

namespace {


static double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Bx - Cx) * (Ay - Cy);
}

/// Point-in-triangle test
/// @note "dot product" test, equivalent to half-plane check
/// @see http://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
static bool point_in_triangle(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy,
        const double& Px, const double& Py ) {

    // Compute signs of dot products
    const bool
            b1 = dot_sign(Px, Py, Ax, Ay, Bx, By) < 0,
            b2 = dot_sign(Px, Py, Bx, By, Cx, Cy) < 0,
            b3 = dot_sign(Px, Py, Cx, Cy, Ax, Ay) < 0;

    // Check if point is in triangle
    // - check for b1 && b2 && b3 only works for triangles ordered counter-clockwise
    // - check for b1==b2 && b2==b3 works for both, equivalent to "all points on the same side"
    return (b1 == b2) && (b2 == b3);
}

}

namespace atlas {
namespace interpolation {
namespace element {

method::Intersect Triag2D::intersects(const Vector2D& p, double edgeEpsilon, double epsilon) const {
    method::Intersect is;
    if( point_in_triangle(v0[0],v0[1],v1[0],v1[1],v2[0],v2[1],p[0],p[1]) ) {
        enum { XX=0, YY=1 };
        const double x0 = v0[XX];
        const double x1 = v1[XX];
        const double x2 = v2[XX];
        const double y0 = v0[YY];
        const double y1 = v1[YY];
        const double y2 = v2[YY];
        const double J = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

        is.u = 1./J * ( (y2 - y0)*p[XX] + (x0 - x2)*p[YY] - x0*y2 + x2*y0 );
        is.v = 1./J * ( (y0 - y1)*p[XX] + (x1 - x0)*p[YY] + x0*y1 - x1*y0 );
        is.success();
    }
    return is;
}

double Triag2D::area() const
{
    NOTIMP;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas

