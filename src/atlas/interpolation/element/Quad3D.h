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

#include <iosfwd>
#include <limits>

#include "atlas/interpolation/Vector3D.h"
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

class Quad3D {
public:
    Quad3D(const double* x0, const double* x1, const double* x2, const double* x3):
        v00(x0), v10(x1), v11(x2), v01(x3) {}

    Quad3D(const PointXYZ& x0, const PointXYZ& x1, const PointXYZ& x2, const PointXYZ& x3):
        Quad3D(x0.data(), x1.data(), x2.data(), x3.data()) {}

    method::Intersect intersects(const method::Ray& r, double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
                                 double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

    bool validate() const;

    double area() const;

    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const Quad3D& p) {
        p.print(s);
        return s;
    }

    const Vector3D& p(int i) {
        if (i == 0)
            return v00;
        if (i == 1)
            return v10;
        if (i == 2)
            return v11;
        if (i == 3)
            return v01;
        throw_OutOfRange("Quad3D::p(i)", i, 4, Here());
    }

private:           // members
    Vector3D v00;  // aka v0
    Vector3D v10;  // aka v1
    Vector3D v11;  // aka v2
    Vector3D v01;  // aka v3
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
