/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/exception/Exceptions.h"

#include "atlas/runtime/Log.h"
#include "atlas/interpolation/element/Quad2D.h"
#include "atlas/interpolation/element/Triag2D.h"

namespace atlas {
namespace interpolation {
namespace element {

//----------------------------------------------------------------------------------------------------------------------

method::Intersect Quad2D::intersects(const Vector2D& p, double edgeEpsilon, double epsilon) const {
    method::Intersect isect; // intersection is false

    Triag2D T013(v00, v10, v01);
    isect = T013.intersects(p, edgeEpsilon, epsilon);
    if(isect)
        return isect;

    Triag2D T231(v11, v01, v10);
    isect = T231.intersects(p, edgeEpsilon, epsilon);
    if(isect) {
        isect.u = 1 - isect.u;
        isect.v = 1 - isect.v;
        return isect;
    }

    return isect.fail();
}

bool Quad2D::validate() const {
    NOTIMP;
}

double Quad2D::area() const {
    Triag2D T013(v00, v10, v01);
    Triag2D T231(v11, v01, v10);

    return T013.area() + T231.area();
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
