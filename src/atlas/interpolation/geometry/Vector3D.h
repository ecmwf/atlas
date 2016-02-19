/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_geometry_Vector3D_h
#define atlas_geometry_Vector3D_h

#include <iostream>
#include <cmath>

#include "atlas/atlas_config.h"

#ifdef ATLAS_HAVE_EIGEN

#define EIGEN_NO_AUTOMATIC_RESIZING
//#define EIGEN_DONT_ALIGN
//#define EIGEN_DONT_VECTORIZE

#include <Eigen/Core>
#include <Eigen/Dense>

#endif

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

#ifdef ATLAS_HAVE_EIGEN

typedef  Eigen::Vector3d Vector3D;

#else

class Vector3D {
  private:

    Vector3D(const double *d): x_(d[0]), y_(d[1]), z_(d[2]) {
    }

    Vector3D(double x, double y, double z): x_(x), y_(y), z_(z) {
    }

  public:

    Vector3D() {
        // Warning, data_ is uninitialised
    }

    static Vector3D Map(const double *data) {
        return Vector3D(data);
    }

    // Operators

    // Vector3D operator*(const Vector3D &) const;
    Vector3D operator-(const Vector3D &other) const {
        return Vector3D(x_ - other.x_, y_ - other.y_, z_ - other.z_);
    }

    Vector3D operator+(const Vector3D &other) const {
        return Vector3D(x_ + other.x_, y_ + other.y_, z_ + other.z_);
    }

    Vector3D operator-() const {
        return Vector3D(-x_, -y_, -z_);
    }

    double norm() const {
        return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    }

    double dot(const Vector3D &other) const {
        return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
    }

    Vector3D cross(const Vector3D &other) const {
        return Vector3D(y_ * other.z_ - z_ * other.y_,
                        z_ * other.x_ - x_ * other.z_,
                        x_ * other.y_ - y_ * other.x_);
    }

    void print(std::ostream &s) const {
        s << "[" << x_ << "," << y_ << "," << z_ << "]";
    }

    friend std::ostream &operator<<(std::ostream &s, const Vector3D &p) {
        p.print(s);
        return s;
    }


  private:
    double x_;
    double y_;
    double z_;

};

Vector3D operator*(double, const Vector3D &);

#endif

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif
