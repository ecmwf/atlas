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

#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/interpolation/method/Intersect.h"
#include "atlas/interpolation/method/Ray.h"

//----------------------------------------------------------------------------------------------------------------------


namespace atlas {
namespace interpolation {
namespace element {

method::Intersect Triag2D::intersects(const method::Ray& r, double edgeEpsilon, double epsilon) const {
    NOTIMP;
}

double Triag2D::area() const
{
    NOTIMP;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas

