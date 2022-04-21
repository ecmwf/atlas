/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iosfwd>
#include <limits>

#include "atlas/interpolation/Vector2D.h"
#include "atlas/interpolation/method/Intersect.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {
struct Ray;
}
namespace element {

//----------------------------------------------------------------------------------------------------------------------

class Quad2D {
public:
    Quad2D(const double* x0, const double* x1, const double* x2, const double* x3):
        v00(x0), v10(x1), v11(x2), v01(x3) {}

    Quad2D(const Point2& x0, const Point2& x1, const Point2& x2, const Point2& x3):
        Quad2D(x0.data(), x1.data(), x2.data(), x3.data()) {}

    Quad2D(const Vector2D& x0, const Vector2D& x1, const Vector2D& x2, const Vector2D& x3):
        Quad2D(x0.data(), x1.data(), x2.data(), x3.data()) {}

    method::Intersect intersects(const Point2& r, double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
                                 double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

    method::Intersect localRemap(const Point2& r, double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
                                 double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

    bool validate() const;

    double area() const;

    const Vector2D& p(int i) {
        if (i == 0) {
            return v00;
        }
        if (i == 1) {
            return v10;
        }
        if (i == 2) {
            return v11;
        }
        if (i == 3) {
            return v01;
        }
        throw_OutOfRange("Quad2D::p(i)", i, 4, Here());
    }

    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const Quad2D& p) {
        p.print(s);
        return s;
    }

private:           // members
    Vector2D v00;  // aka v0
    Vector2D v10;  // aka v1
    Vector2D v11;  // aka v2
    Vector2D v01;  // aka v3

    bool inQuadrilateral(const Vector2D& p, double tolerance = 5 * std::numeric_limits<double>::epsilon()) const;

    static double cross2d(const Vector2D& a, const Vector2D& b) { return a.x() * b.y() - a.y() * b.x(); }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
