/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Vector2D.h"

#if !ATLAS_HAVE_EIGEN

#include <iostream>

namespace atlas {
namespace interpolation {

void Vector2D::print(std::ostream& s) const {
    s << "[" << x() << "," << y() << "]";
}

}  // namespace interpolation
}  // namespace atlas
#endif
