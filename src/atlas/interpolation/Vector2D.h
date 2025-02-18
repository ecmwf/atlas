/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/config.h"

#if ATLAS_HAVE_EIGEN

#define EIGEN_NO_AUTOMATIC_RESIZING
//#define EIGEN_DONT_ALIGN
//#define EIGEN_DONT_VECTORIZE

ATLAS_SUPPRESS_WARNINGS_PUSH
ATLAS_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE
ATLAS_SUPPRESS_WARNINGS_CODE_IS_UNREACHABLE
#include <Eigen/Core>
#include <Eigen/Dense>
ATLAS_SUPPRESS_WARNINGS_POP
#else

#include <cmath>
#include <iosfwd>

#endif

namespace atlas {
namespace interpolation {

//----------------------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_EIGEN

typedef Eigen::Vector2d Vector2D;

#else

class Vector2D {
public:
    Vector2D(const double* d) {
        xy_[0] = d[0];
        xy_[1] = d[1];
    }

    Vector2D(double x, double y) {
        xy_[0] = x;
        xy_[1] = y;
    }

    Vector2D() {
        // Warning, data_ is uninitialised
    }

    static Vector2D Map(const double* data) { return Vector2D(data); }

    // Operators

    double operator[](size_t i) const { return xy_[i]; }

    // Vector2D operator*(const Vector2D &) const;
    Vector2D operator-(const Vector2D& other) const { return Vector2D(x() - other.x(), y() - other.y()); }

    Vector2D operator+(const Vector2D& other) const { return Vector2D(x() + other.x(), y() + other.y()); }

    Vector2D operator-() const { return Vector2D(-x(), -y()); }

    double norm() const { return sqrt(squaredNorm()); }

    double squaredNorm() const { return x() * x() + y() * y(); }

    double dot(const Vector2D& other) const { return x() * other.x() + y() * other.y(); }

    double cross(const Vector2D& other) const { return x() * other.y() - y() * other.x(); }

    void print(std::ostream& s) const;

    friend std::ostream& operator<<(std::ostream& s, const Vector2D& p) {
        p.print(s);
        return s;
    }

    double* data() { return xy_; }

    const double* data() const { return xy_; }

    double x() const { return xy_[0]; }
    double y() const { return xy_[1]; }

private:
    double xy_[2];
};

Vector2D operator*(double, const Vector2D&);

#endif

//----------------------------------------------------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace atlas
