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
#include <utility>

#include "atlas/util/Constants.h"
#include "atlas/util/Point.h"

#include "eckit/geometry/UnitSphere.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

using eckit::geometry::UnitSphere;

/// @brief   Calculate great-cricle course between points
///
/// @details Calculates the direction (clockwise from north) of a great-circle
///          arc between lonLat1 and lonLat2. Returns the direction of the arc
///          at lonLat1 (first) and lonLat2 (second). Angle is normalised to the
///          range of atan2 (usually (-180, 180]). All input and output values
///          are in units of degrees.
/// @ref     https://en.wikipedia.org/wiki/Great-circle_navigation
///
inline std::pair<double, double> greatCircleCourse(const Point2& lonLat1,
                                                   const Point2& lonLat2) {

  const auto lambda1 = lonLat1[0] * Constants::degreesToRadians();
  const auto lambda2 = lonLat2[0] * Constants::degreesToRadians();
  const auto phi1 = lonLat1[1] * Constants::degreesToRadians();
  const auto phi2 = lonLat2[1] * Constants::degreesToRadians();

  const auto sinLambda12 = std::sin(lambda2 - lambda1);
  const auto cosLambda12 = std::cos(lambda2 - lambda1);
  const auto sinPhi1 = std::sin(phi1);
  const auto sinPhi2 = std::sin(phi2);
  const auto cosPhi1 = std::cos(phi1);
  const auto cosPhi2 = std::cos(phi2);

  const auto alpha1 =
      std::atan2(cosPhi2 * sinLambda12,
                 cosPhi1 * sinPhi2 - sinPhi1 * cosPhi2 * cosLambda12);

  const auto alpha2 =
      std::atan2(cosPhi1 * sinLambda12,
                 -cosPhi2 * sinPhi1 + sinPhi2 * cosPhi1 * cosLambda12);

  return std::make_pair(alpha1 * Constants::radiansToDegrees(),
                        alpha2 * Constants::radiansToDegrees());
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
