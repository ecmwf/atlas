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
///           "Benchmarking Regridding Libraries Used in Earth System Modelling"
///           by Sophie Valcke, Andreas Piacentini, Gabriel Jonville, MDPI 2022
///           as the sinusoid analytical function in Sec 2.1.2.
///           The longitude (lon) and latitude (lat) are assumed to be in degrees,
///
double MDPI_sinusoid(double lon, double lat);
double MDPI_harmonic(double lon, double lat);
double MDPI_vortex(double lon, double lat);
double MDPI_gulfstream(double lon, double lat);

extern "C" {
    double atlas__functions__MDPI_sinusoid(double& lon, double& lat);
    double atlas__functions__MDPI_harmonic(double& lon, double& lat);
    double atlas__functions__MDPI_vortex(double& lon, double& lat);
    double atlas__functions__MDPI_gulfstream(double& lon, double& lat);
}

}  // namespace function

}  // namespace util

}  // namespace atlas
