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

/// \brief An analytic function that provides a vortex rollup on a 2D sphere
///
/// \detailed The formula is found in
///           "A Lagrangian Particle Method with Remeshing for Tracer Transport on the Sphere"
///           by Peter Bosler, James Kent, Robert Krasny, Christiane Jablonowski, JCP 2015
///           as the tracer density in Test Case 1.
///           The longitude (lon) and latitude (lat) are assumed to be in degrees,
///           The time parameter associated with the vortex rollup is set by (t).
///
///           The time it takes for the counter-rotating vortices along
///           the equator to return to its original position takes
///           time t = 1.0;
///
double vortex_rollup(double lon, double lat, double t);

}  // namespace function

}  // namespace util

}  // namespace atlas
