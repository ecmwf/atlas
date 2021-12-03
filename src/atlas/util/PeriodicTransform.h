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

#include <array>
#include "atlas/util/LonLatMicroDeg.h"

namespace atlas {
namespace util {

class PeriodicTransform {
protected:
    double x_translation_;

public:
    PeriodicTransform() { x_translation_ = 360.; }
    PeriodicTransform(const double& translation) { x_translation_ = translation; }

    void operator()(double source[2], double dest[2], double direction, double scale = 1.) const {
        dest[0] = source[0] + direction * x_translation_ * scale;
        dest[1] = source[1];
    }

    void operator()(int source[2], int dest[2], int direction, int scale = 1) const {
        dest[0] = source[0] + direction * static_cast<int>(x_translation_ * scale);
        dest[1] = source[1];
    }

    void operator()(double inplace[2], double direction, double scale = 1.) const {
        inplace[0] = inplace[0] + direction * x_translation_ * scale;
        // inplace[1] = inplace[1]; null operation
    }

    void operator()(int inplace[2], int direction, int scale = 1) const {
        inplace[0] = inplace[0] + direction * static_cast<int>(x_translation_ * scale);
        // inplace[1] = inplace[1]; null operation
    }

    void operator()(util::LonLatMicroDeg& inplace, int direction) const {
        inplace.set_lon(inplace.lon() + direction * util::microdeg(x_translation_));
        // inplace.set_lat( inplace.lat() ); null operation
    }

    void operator()(std::array<double, 2>& inplace) const {
        inplace[0] = inplace[0] + x_translation_;
        // inplace[1] = inplace[1]; null operation
    }
};

}  // namespace util
}  // namespace atlas
