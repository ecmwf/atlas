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

typedef Eigen::Vector3d Vector3D;

#else

class Vector3D {
public:
    Vector3D(const double* d) {
        xyz_[0] = d[0];
        xyz_[1] = d[1];
        xyz_[2] = d[2];
    }

    Vector3D(double x, double y, double z) {
        xyz_[0] = x;
        xyz_[1] = y;
        xyz_[2] = z;
    }

    Vector3D() {
        // Warning, data_ is uninitialised
    }

    static Vector3D Map(const double* data) { return Vector3D(data); }

    // Operators

    double operator[](size_t i) const { return xyz_[i]; }

    // Vector3D operator*(const Vector3D &) const;
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x() - other.x(), y() - other.y(), z() - other.z());
    }

    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x() + other.x(), y() + other.y(), z() + other.z());
    }

    Vector3D operator-() const { return Vector3D(-x(), -y(), -z()); }

    Vector3D operator/(double a) const { return Vector3D(x() / a, y() / a, z() / a); }

    Vector3D operator*(double a) const { return Vector3D(x() * a, y() * a, z() * a); }

    double norm() const { return std::sqrt(squaredNorm()); }

    double squaredNorm() const { return x() * x() + y() * y() + z() * z(); }

    double dot(const Vector3D& other) const { return x() * other.x() + y() * other.y() + z() * other.z(); }

    Vector3D cross(const Vector3D& other) const {
        return Vector3D(y() * other.z() - z() * other.y(), z() * other.x() - x() * other.z(),
                        x() * other.y() - y() * other.x());
    }

    void print(std::ostream& s) const;

    friend std::ostream& operator<<(std::ostream& s, const Vector3D& p) {
        p.print(s);
        return s;
    }

    double* data() { return xyz_; }

    const double* data() const { return xyz_; }

private:
    double x() const { return xyz_[0]; }
    double y() const { return xyz_[1]; }
    double z() const { return xyz_[2]; }
    double xyz_[3];
};

Vector3D operator*(double, const Vector3D&);

#endif

//----------------------------------------------------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace atlas
