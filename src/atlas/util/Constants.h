/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_Constants_h
#define atlas_util_Constants_h

#include <cmath>
#include "atlas/internals/atlas_defines.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

/// Some usefull constants
/// @note These could be static const constants, but then the initialization would be in the .cc
///       which would preclude constant optimization.
///       With C++11 constexpr the constants can be initialized in the class header.
struct Constants
{
    static constexpr double radiansToDegrees() { return 180. * M_1_PI; }
    static constexpr double degreesToRadians() { return M_PI / 180.; }

};

//------------------------------------------------------------------------------------------------------

struct Earth
{
    static constexpr double radiusInMeters() { return 6371229.; }
    static constexpr double radiusInKm()     { return radiusInMeters() / 1.0e3; }

    static constexpr double areaInSqMeters() { return 4. * M_PI * radiusInMeters() * radiusInMeters(); }
    static constexpr double areaInSqKm()     { return 4. * M_PI * radiusInKm()     * radiusInKm();     }
};

//------------------------------------------------------------------------------------------------------

} // namespace util
} // namespace atlas

#endif // atlas_util_Constants_h
