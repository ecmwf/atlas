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
#include <iosfwd>

#include "atlas/util/Point.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace util {

// Compute coordinates of a point on a rotated sphere with specified pole
class Rotation {
public:
    Rotation(const PointLonLat& south_pole, double rotation_angle = 0.);
    Rotation(const eckit::Parametrisation&);

    bool rotated() const { return rotated_; }

    PointLonLat southPole() const { return spole_; }
    PointLonLat northPole() const { return npole_; }
    double rotationAngle() const { return angle_; }

    PointLonLat rotate(const PointLonLat&) const;
    PointLonLat unrotate(const PointLonLat&) const;

    void rotate(double crd[]) const;
    void unrotate(double crd[]) const;

private:
    void precompute();

    void print(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Rotation&);

private:
    PointLonLat npole_ = {-180., 90.};  // North Pole
    PointLonLat spole_ = {0., -90.};    // South Pole
    double angle_      = {0.};

    using RotationMatrix = std::array<std::array<double, 3>, 3>;

    RotationMatrix rotate_;    // rotate   matrix
    RotationMatrix unrotate_;  // unrotate matrix

    bool rotation_angle_only_;
    bool rotated_;
};

}  // namespace util
}  // namespace atlas
