/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "atlas/interpolation/Vector2D.h"
#include "atlas/interpolation/element/Quad2D.h"
#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace interpolation {
namespace element {

//----------------------------------------------------------------------------------------------------------------------

method::Intersect Quad2D::intersects(const PointXY& r, double edgeEpsilon, double epsilon) const {
    method::Intersect isect;  // intersection is false

    /* Split quadrilateral into two triangles, points are labelled counter-clockwise.
     *   v01-----------v11
     *    /  " .        `.
     *   /         " .    `.
     *  /                " .`.
     * v00-------------------v10
     *
     */

    Triag2D T013(v00, v10, v01);
    isect = T013.intersects(r, edgeEpsilon, epsilon);
    if (isect) {
        return isect;
    }

    Triag2D T231(v11, v01, v10);
    isect = T231.intersects(r, edgeEpsilon, epsilon);
    if (isect) {
        isect.u = 1 - isect.u;
        isect.v = 1 - isect.v;
        return isect;
    }

    return isect.fail();
}

method::Intersect Quad2D::localRemap(const PointXY& p, double edgeEpsilon, double epsilon) const {
    method::Intersect isect;

    // work out if point is within the polygon
    if (!inQuadrilateral({p.x(), p.y()}))
        return isect.fail();

    // calculate local coordinates via remapping technique
    // if that fails fall back to finite-element scheme.
    //
    // ax**2 + bx + x = 0
    auto solve_quadratic = [](const double a, const double b, const double c) {
        double det = b * b - 4. * a * c;
        Roots roots;
        if (det >= 0.) {
            double inv_two_a = 1. / (2. * a);
            double sqrt_det = std::sqrt(det);
            roots.a = (-b + sqrt_det) * inv_two_a;
            roots.b = (-b - sqrt_det) * inv_two_a;
        }
        return roots;
    };

    // ax + b = 0
    auto solve_linear = [](const double a, const double b) { return -b / a; };

    auto validWeight = [&](const double w) { return (w > - epsilon) && (w < (1. + epsilon)); };

    // solve for u and v where:
    // w1 = ( 1 - u ) * ( 1 - v )
    // w2 = u * ( 1 - v )
    // w3 = u * v
    // w4 = ( 1 - u ) * v

    Vector2D ray(p.x(), p.y());
    Vector2D vA = v00 - ray;
    Vector2D vB = v10 - v00;
    Vector2D vC = v01 - v00;
    Vector2D vD = v00 - v10 - v01 + v11;

    // solve for v
    double a = cross2d(vC, vD);
    double b = cross2d(vC, vB) + cross2d(vA, vD);
    double c = cross2d(vA, vB);

    if (abs(a) > epsilon) {
        Roots roots = solve_quadratic(a, b, c);
        if (validWeight(roots.a)) {
            isect.v = roots.a;
        }
        else if (validWeight(roots.b)) {
            isect.v = roots.b;
        }
        else {
            isect.v = BAD_WEIGHT_VALUE;
        }
    }
    else if (abs(b) > epsilon) {
        isect.v = solve_linear(b, c);
    }
    else {
        isect.v = BAD_WEIGHT_VALUE;
    }

    // solve for u
    a = cross2d(vB, vD);
    b = cross2d(vB, vC) + cross2d(vA, vD);
    c = cross2d(vA, vC);

    if (abs(a) > epsilon) {
        Roots roots = solve_quadratic(a, b, c);
        if (validWeight(roots.a)) {
            isect.u = roots.a;
        }
        else if (validWeight(roots.b)) {
            isect.u = roots.b;
        }
        else {
            isect.u = BAD_WEIGHT_VALUE;
        }
    }
    else if (abs(b) > epsilon) {
        isect.u = solve_linear(b, c);
    }
    else {
        isect.u = BAD_WEIGHT_VALUE;
    }

    if ((isect.u == BAD_WEIGHT_VALUE) || (isect.v == BAD_WEIGHT_VALUE)) {
        return isect.fail();
    }
    else {
        return isect.success();
    }
}

bool Quad2D::validate() const {
    // normal for sub-triangle T231

    Vector2D E23 = v01 - v11;
    Vector2D E21 = v10 - v11;

    double N231 = cross2d(E23, E21);

    // normal for sub-triangle T013

    Vector2D E01 = v10 - v00;
    Vector2D E03 = v01 - v00;

    double N013 = cross2d(E01, E03);

    // normal for sub-triangle T120

    Vector2D E12 = -E21;
    Vector2D E10 = -E01;

    double N120 = cross2d(E12, E10);

    // normal for sub-triangle T302

    Vector2D E30 = -E03;
    Vector2D E32 = -E23;

    double N302 = cross2d(E30, E32);

    // all normals must point same way

    double dot02 = N231 * N013;
    double dot23 = N013 * N120;
    double dot31 = N120 * N302;
    double dot10 = N302 * N231;

    // all normals must point same way
    bool is_inside = ((dot02 >= 0. && dot23 >= 0. && dot31 >= 0. && dot10 >= 0.) ||
                      (dot02 <= 0. && dot23 <= 0. && dot31 <= 0. && dot10 <= 0.));
    return is_inside;
}

double Quad2D::area() const {
    return std::abs(0.5 * cross2d((v01 - v00), (v11 - v00))) + std::abs(0.5 * cross2d((v10 - v11), (v01 - v11)));
}

bool Quad2D::inQuadrilateral(const Vector2D& p) const {
    // point p must be on the inside of all quad edges to be inside the quad.
    double zst1 = cross2d(p - v00, p - v10);
    if (zst1 >= 0.0) {
        double zst2 = cross2d(p - v10, p - v11);
        if (zst2 >= 0.0) {
            double zst3 = cross2d(p - v11, p - v01);
            if (zst3 >= 0.0) {
                double zst4 = cross2d(p - v01, p - v00);
                if (zst4 >= 0.0)
                    return true;
            }
        }
    }

    return false;
}

void Quad2D::print(std::ostream& s) const {
    auto printVector2D = [&s](const Vector2D& v) { s << "[" << v[0] << "," << v[1] << "]"; };
    s << "Quad2D[";
    printVector2D(v00);
    s << ", ";
    printVector2D(v10);
    s << ", ";
    printVector2D(v11);
    s << ", ";
    printVector2D(v01);
    s << "]";
}


//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
