/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/interpolation/method/Intersect.h"

#include <iostream>

namespace atlas {
namespace interpolation {
namespace method {

Intersect::Intersect(): u(0.), v(0.), t(0.), success_(false) {}

void Intersect::print(std::ostream& s) const {
    s << "Intersect[u=" << u << ",v=" << v << ",t=" << t << ",success=" << success_ << "]";
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
