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

#include "atlas/interpolation/Vector3D.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

/// Ray trace data structure

struct Ray {
    Vector3D orig;
    Vector3D dir;

    /// initializes ray with origin in point and direction to (0,0,0)
    explicit Ray(const double* p);

    Ray(const double* o, const double* d);

    explicit Ray(const PointXYZ&);

    Ray(const PointXYZ&, const Vector3D&);

    Vector3D operator()(double t) const { return orig + t * dir; }

    void print(std::ostream& s) const;

    friend std::ostream& operator<<(std::ostream& s, const Ray& p);
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
