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

#include <cassert>

namespace atlas {
namespace util {

inline int microdeg(const double& deg) {
    assert(deg < 2145.);   // Since INT_MAX ==  2147483647
    assert(deg > -2145.);  // Since INT_MIN == –2147483648
    return static_cast<int>(deg < 0 ? deg * 1.e6 - 0.5 : deg * 1.e6 + 0.5);
}

}  // namespace util
}  // namespace atlas
