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

#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <string>

#include "eckit/geometry/Point2.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {

using Point2 = eckit::geometry::Point2;

//---------------------------------------------------------------------------------------------------------------------
namespace projection {

class Jacobian : public std::array<std::array<double, 2>, 2> {
public:
    using std::array<std::array<double, 2>, 2>::array;

    static Jacobian identity() { return Jacobian{1., 0., 0., 1.}; }

    Jacobian() = default;

    Jacobian(double j00, double j01, double j10, double j11) {
        (*this)[0][0] = j00;
        (*this)[0][1] = j01;
        (*this)[1][0] = j10;
        (*this)[1][1] = j11;
    }

    Jacobian(std::initializer_list<std::initializer_list<double>> list):
        Jacobian{*(list.begin()->begin()), *(list.begin()->begin() + 1), *((list.begin() + 1)->begin()),
                 *((list.begin() + 1)->begin() + 1)} {}

    Jacobian operator-(const Jacobian& jac) const {
        return Jacobian{(*this)[0][0] - jac[0][0], (*this)[0][1] - jac[0][1], (*this)[1][0] - jac[1][0],
                        (*this)[1][1] - jac[1][1]};
    }

    Jacobian operator+(const Jacobian& jac) const {
        return Jacobian{(*this)[0][0] + jac[0][0], (*this)[0][1] + jac[0][1], (*this)[1][0] + jac[1][0],
                        (*this)[1][1] + jac[1][1]};
    }

    Jacobian operator*(double a) const {
        return Jacobian{(*this)[0][0] * a, (*this)[0][1] * a, (*this)[1][0] * a, (*this)[1][1] * a};
    }

    Jacobian operator*(const Jacobian& jac) const {
        return Jacobian{(*this)[0][0] * jac[0][0] + (*this)[0][1] * jac[1][0],
                        (*this)[0][0] * jac[0][1] + (*this)[0][1] * jac[1][1],
                        (*this)[1][0] * jac[0][0] + (*this)[1][1] * jac[1][0],
                        (*this)[1][0] * jac[0][1] + (*this)[1][1] * jac[1][1]};
    }

    Point2 operator*(const Point2& x) const {
        return Point2{(*this)[0][0] * x[0] + (*this)[0][1] * x[1], (*this)[1][0] * x[0] + (*this)[1][1] * x[1]};
    }

    double norm() const {
        return std::sqrt((*this)[0][0] * (*this)[0][0] + (*this)[0][1] * (*this)[0][1] + (*this)[1][0] * (*this)[1][0] +
                         (*this)[1][1] * (*this)[1][1]);
    }

    double determinant() const { return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0]; }

    Jacobian inverse() const {
        return Jacobian{(*this)[1][1], -(*this)[0][1], -(*this)[1][0], (*this)[0][0]} * (1. / determinant());
    }

    Jacobian transpose() const { return Jacobian{(*this)[0][0], (*this)[1][0], (*this)[0][1], (*this)[1][1]}; }

    double dx_dlon() const { return (*this)[JDX][JDLON]; }
    double dy_dlon() const { return (*this)[JDY][JDLON]; }
    double dx_dlat() const { return (*this)[JDX][JDLAT]; }
    double dy_dlat() const { return (*this)[JDY][JDLAT]; }

    double dlon_dx() const { return (*this)[JDLON][JDX]; }
    double dlon_dy() const { return (*this)[JDLON][JDY]; }
    double dlat_dx() const { return (*this)[JDLAT][JDX]; }
    double dlat_dy() const { return (*this)[JDLAT][JDY]; }

    friend std::ostream& operator<<(std::ostream& os, const Jacobian& jac) {
        os << jac[0][0] << " " << jac[0][1] << "\n" << jac[1][0] << " " << jac[1][1];
        return os;
    }

private:
    enum
    {
        JDX = 0,
        JDY = 1
    };
    enum
    {
        JDLON = 0,
        JDLAT = 1
    };
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace projection
}  // namespace atlas
