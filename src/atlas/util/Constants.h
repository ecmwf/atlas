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

#include <cmath>

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

/// Some useful constants
struct Constants {
    static constexpr double radiansToDegrees() { return 180. * M_1_PI; }
    static constexpr double degreesToRadians() { return M_PI / 180.; }
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
