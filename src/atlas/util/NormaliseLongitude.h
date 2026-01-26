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

/// @file NormaliseLongitude.h

namespace atlas {
namespace util {

class NormaliseLongitude {
public:
    // Normalise longitude between (west - eps, east - eps ) with west = 0., east = 360.
    constexpr NormaliseLongitude(): west_(-eps_), east_(360. - eps_) {}

    // Normalise longitude between ( west-eps, east-eps )  with east = west + 360
    constexpr NormaliseLongitude(double west): west_(west - eps_), east_(west + 360. - eps_) {}

    // Normalise longitude between ( west-eps, east+eps )
    constexpr NormaliseLongitude(double west, double east): west_(west - eps_), east_(east + eps_) {}

    double operator()(double lon) const {
        while (lon < west_) {
            lon += 360.;
        }
        while (lon > east_) {
            lon -= 360.;
        }
        return lon;
    }

private:
    double west_;
    double east_;

public:
    static constexpr double eps_ = 1.e-11;
};

}  // namespace util
}  // namespace atlas
