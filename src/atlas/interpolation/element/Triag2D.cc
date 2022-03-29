/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "atlas/interpolation/Vector2D.h"
#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace interpolation {
namespace element {

//----------------------------------------------------------------------------------------------------------------------

method::Intersect Triag2D::intersects(const PointXY& r, double edgeEpsilon, double epsilon) const {
    method::Intersect isect;  // intersection is false

    /*
     *                ,v11
     *     e2     , ^   `.
     *        , ^     r   `.
     *    , ^               `.
     * v00-------------------v10
     *             e1
     */

    Vector2D rvec{r.data()};

    if (!inTriangle(rvec, epsilon)) {
        return isect.fail();
    }

    Vector2D e1{v10 - v00};
    Vector2D e2{v11 - v00};
    Vector2D pvec{rvec - v00};

    // solve u e1 + v e2 = pvec for u and v
    float invDet = 1. / (e1.x() * e2.y() - e2.x() * e1.y());
    isect.u      = (pvec.x() * e2.y() - e2.x() * pvec.y()) * invDet;
    isect.v      = (e1.x() * pvec.y() - pvec.x() * e1.y()) * invDet;

    // clamp values between 0 and 1
    if (isect.u > -edgeEpsilon || isect.u < 1. + edgeEpsilon ||
        isect.v > -edgeEpsilon || isect.v < 1. + edgeEpsilon) {

        return isect.success();
    }
    else {
        return isect.fail();
    }
}

bool Triag2D::validate() const {
    // normal for sub-triangle T012

    Vector2D E01 = v10 - v00;
    Vector2D E02 = v11 - v00;

    double N012 = cross2d(E01, E02);

    // normal for sub-triangle T120

    Vector2D E12 = v11 - v10;
    Vector2D E10 = v00 - v10;

    double N120 = cross2d(E12, E10);

    // normal for sub-triangle T201

    Vector2D E20 = -E02;
    Vector2D E21 = -E12;

    double N201 = cross2d(E20, E21);

    // all normals must point same way

    double dot1 = N120 * N012;
    double dot2 = N201 * N120;
    double dot3 = N201 * N012;

    // all normals must point same way
    bool is_inside = ((dot1 >= 0. && dot2 >= 0. && dot3 >= 0.) || (dot1 <= 0. && dot2 <= 0. && dot3 <= 0.));
    return is_inside;
}

double Triag2D::area() const {
    return std::abs(0.5 * cross2d((v10 - v00), (v11 - v00)));
}

bool Triag2D::inTriangle(const Vector2D& p, double epsilon) const {
    // point p must be on the inside of all triangle edges to be inside the triangle.
    return cross2d(p - v00, p - v10) > -epsilon &&
           cross2d(p - v10, p - v11) > -epsilon &&
           cross2d(p - v11, p - v00) > -epsilon;
}

void Triag2D::print(std::ostream& s) const {
    auto printVector2D = [&s](const Vector2D& v) { s << "[" << v[0] << "," << v[1] << "]"; };
    s << "Triag2D[";
    printVector2D(v00);
    s << ", ";
    printVector2D(v11);
    s << ", ";
    printVector2D(v10);
    s << "]";
}


//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
