/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iostream>

#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

#include "atlas/util/function/SlottedCylinder_function.h"

namespace atlas {

namespace util {

namespace function {

double SlottedCylinder(double lon, double lat, double scale) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();
    double x = lon - M_PI;
    double y = lat;
    double r2 = x * x + y * y;
    if (r2 <= 1.5 * scale && (std::abs(x) >= 0.25 * scale || y >= scale)) {
        return 1.;
    }
    return 0.;
}


extern "C" {
    double atlas__functions__SlottedCylinder(double& lon, double& lat, double& scale) {
        return SlottedCylinder(lon, lat, scale);
    }
}

}  // namespace function
}  // namespace util
}  // namespace atlas
