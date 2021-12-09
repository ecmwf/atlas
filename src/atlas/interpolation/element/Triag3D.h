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

/// Triangle structure
/// Implements @link
/// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm

class Triag3D {
public:  // types
    Triag3D(const double* x0, const double* x1, const double* x2): v0(x0), v1(x1), v2(x2) {}

    Triag3D(const PointXYZ& x0, const PointXYZ& x1, const PointXYZ& x2): Triag3D(x0.data(), x1.data(), x2.data()) {}

    Triag3D(const Vector3D& x0, const Vector3D& x1, const Vector3D& x2): Triag3D(x0.data(), x1.data(), x2.data()) {}

    method::Intersect intersects(const method::Ray& r, double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
                                 double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

    double area() const;

    void print(std::ostream& s) const;

    friend std::ostream& operator<<(std::ostream& s, const Triag3D& p) {
        p.print(s);
        return s;
    }

    const Vector3D& p(int i) {
        if (i == 0)
            return v0;
        if (i == 1)
            return v1;
        if (i == 2)
            return v2;
        throw_OutOfRange("Triag3D::p(i)", i, 3, Here());
    }

private:  // members
    Vector3D v0;
    Vector3D v1;
    Vector3D v2;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
