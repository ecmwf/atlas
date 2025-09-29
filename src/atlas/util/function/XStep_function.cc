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

#include "atlas/util/function/MDPI_functions.h"

namespace atlas {

namespace util {

namespace function {

double XStep(double lon, double lat) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();
    double x1 = lon - M_PI;
    double x2 = -lon + M_PI;
    if ((lat > x1) && (lat < x2)) {
        return 1.;
    }
    if ((lat < x1) && (lat > x2)) {
        return -1.;
    }
    return 0.;
}


extern "C" {
    double atlas__functions__XStep(double& lon, double& lat) {
        return XStep(lon, lat);
    }
}

}  // namespace function
}  // namespace util
}  // namespace atlas
