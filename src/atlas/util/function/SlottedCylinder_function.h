/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

namespace atlas {

namespace util {

namespace function {

/// \brief Analytic functions from MDPI paper on regridding
///
/// \detailed The formula is found in
///           "Fully Multidimensional Flux-Corrected Transport Algorthim for Fluids"
///           by Steven T. Zalesak, JCP 1979
///           as given in Fig. 11
///           The longitude (lon) and latitude (lat) are assumed to be in radians.
///
double SlottedCylinder(double lon, double lat);


extern "C" {
    double atlas__functions__SlottedCylinder(double& lon, double& lat);
}

}  // namespace function

}  // namespace util

}  // namespace atlas
