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
namespace element {

static constexpr double BAD_WEIGHT_VALUE = -1.;

//----------------------------------------------------------------------------------------------------------------------

class Triag2D {
public:
    Triag2D(const double* x0, const double* x1, const double* x2): v00(x0), v10(x1), v11(x2) {}

    Triag2D(const PointXY& x0, const PointXY& x1, const PointXY& x2): Triag2D(x0.data(), x1.data(), x2.data()) {}

    Triag2D(const Vector2D& x0, const Vector2D& x1, const Vector2D& x2): Triag2D(x0.data(), x1.data(), x2.data()) {}

    method::Intersect intersects(const PointXY& r, double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
                                 double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

    bool validate() const;

    double area() const;

    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const Triag2D& p) {
        p.print(s);
        return s;
    }

    const Vector2D& p(int i) {
        if (i == 0)
            return v00;
        if (i == 1)
            return v10;
        if (i == 2)
            return v11;
        throw_OutOfRange("Triag2D::p(i)", i, 3, Here());
    }

private:           // members
    Vector2D v00;  // aka v0
    Vector2D v10;  // aka v1
    Vector2D v11;  // aka v2

    static double cross2d(const Vector2D& a, const Vector2D& b) { return a.x() * b.y() - a.y() * b.x(); }

    bool inTriangle(const Vector2D& p, double tolerance = 5 * std::numeric_limits<double>::epsilon()) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
