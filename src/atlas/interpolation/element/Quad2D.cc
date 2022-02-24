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

namespace  {

constexpr double badRoot = std::numeric_limits<double>::lowest();

struct Roots {
    double x1;
    double x2;
    double epsilon;
};


// Solve ax^2 + bx + c = 0 for x
Roots solveQuadratic(double a, double b, double c, double epsilon, double rootEpsilon) {

    if (std::abs(a) < epsilon ) {
        const auto x = -c / b;
        return Roots{x, x, rootEpsilon};
    }
    else {

        // x = (-b +/- sqrt(b^2 - 4ac)) / (2a)

        const auto det = b * b - 4. * a * c;
        if (det >= -epsilon) {

            const auto sqrtDet = std::sqrt(std::max(0., det));

            const auto inv2a = 0.5 / a;

            const auto f1 = -b * inv2a;
            const auto f2 = sqrtDet * inv2a;

            // Scale root epsilon by larger term.
            rootEpsilon *= std::max({1., std::abs(f1), std::abs(f2)});

            return Roots{f1 + f2, f1 - f2, rootEpsilon};

        }
        else {
            return Roots{badRoot, badRoot, rootEpsilon};
        }

    }
}

bool validRoot(double root, double rootEpsilon) {
    return root > -rootEpsilon && root < 1. + rootEpsilon;
}

}

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

    // Perform "inverse" bilinear interpolation.
    // Method described by Inigo Quilez:
    // https://www.iquilezles.org/www/articles/ibilinear/ibilinear.htm
    // This implementation follows the notation given on wikipedia (01/02/2022):
    // https://en.wikipedia.org/wiki/Bilinear_interpolation#Inverse_and_generalization

    method::Intersect isect{};

    // work out if point is within the polygon
    const Vector2D point = Vector2D(p.data());
    if (!inQuadrilateral(point, epsilon)) {
        return isect.fail();
    }

    // Get helper vectors.
    const Vector2D A = v00 - point;
    const Vector2D B = v10 - v00;
    const Vector2D C = v01 - v00;
    const Vector2D D = v00 - v10 + v11 - v01;

    // Get helper scalars.
    const double a = cross2d(A, B);
    const double b = cross2d(A, C);
    const double c = cross2d(A, D);
    const double d = cross2d(B, C);
    const double e = cross2d(B, D);
    const double f = cross2d(C, D);


    // Solve for u.
    const auto uRoots = solveQuadratic(e, c + d, b, epsilon, edgeEpsilon);
    isect.u = validRoot(uRoots.x1, uRoots.epsilon) ? uRoots.x1 : uRoots.x2;

    // Solve for v.
    const auto vRoots = solveQuadratic(f, c - d, a, epsilon, edgeEpsilon);
    isect.v = validRoot(vRoots.x1, vRoots.epsilon) ? vRoots.x1 : vRoots.x2;

    return validRoot(isect.u, uRoots.epsilon) && validRoot(isect.u, vRoots.epsilon) ?
                isect.success() : isect.fail();
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
    bool is_inside = (dot02 >= 0. && dot23 >= 0. && dot31 >= 0. && dot10 >= 0.) ||
                     (dot02 <= 0. && dot23 <= 0. && dot31 <= 0. && dot10 <= 0.);
    return is_inside;
}

double Quad2D::area() const {
    return std::abs(0.5 * cross2d((v01 - v00), (v11 - v00))) + std::abs(0.5 * cross2d((v10 - v11), (v01 - v11)));
}

bool Quad2D::inQuadrilateral(const Vector2D& p, double epsilon) const {
    // point p must be on the inside of all quad edges to be inside the quad.
    const double tol = -epsilon;
    return cross2d(p - v00, p - v10) > tol &&
           cross2d(p - v10, p - v11) > tol &&
           cross2d(p - v11, p - v01) > tol &&
           cross2d(p - v01, p - v00) > tol;
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
