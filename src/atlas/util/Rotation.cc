/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Rotation.h"

#include <cmath>
#include <iostream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/UnitSphere.h"

namespace atlas {
namespace util {

namespace {

PointLonLat wrap_latitude(const PointLonLat& p) {
    double lon = p.lon();
    double lat = p.lat();

    while (lat > 90.) {
        lon += 180.;
        lat = 180. - lat;
    }

    while (lat < -90.) {
        lon -= 180.;
        lat = -180. - lat;
    }

    return {lon, lat};
}

double wrap_angle(double a) {
    PointLonLat angle{a, 0};
    angle.normalise();
    return angle.lon();
}

}  // namespace

void Rotation::print(std::ostream& out) const {
    out << "north_pole:" << npole_ << ", south_pole:" << spole_ << ", rotation_angle:" << angle_;
}
std::ostream& operator<<(std::ostream& out, const Rotation& r) {
    r.print(out);
    return out;
}

PointLonLat Rotation::rotate(const PointLonLat& p) const {
    PointLonLat rotated(p);
    rotate(rotated.data());
    return rotated;
}

PointLonLat Rotation::unrotate(const PointLonLat& p) const {
    PointLonLat unrotated(p);
    unrotate(unrotated.data());
    return unrotated;
}

void Rotation::precompute() {
    const double theta = -(90. + spole_.lat()) * Constants::degreesToRadians();
    const double phi   = -spole_.lon() * Constants::degreesToRadians();

    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    const double sin_phi   = std::sin(phi);
    const double cos_phi   = std::cos(phi);

    auto eq = [](double a, double b) { return eckit::types::is_approximately_equal(a, b, 1.e-12); };

    rotated_ = !eq(angle_, 0) || !eq(sin_theta, 0) || !eq(cos_theta, 1) || !eq(sin_phi, 0) || !eq(cos_phi, 1);
    rotation_angle_only_ = eq(sin_theta, 0) && eq(cos_theta, 1) && rotated_;

    if (!rotated_) {
        rotate_ = unrotate_ = RotationMatrix{{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}};
        return;
    }

    // Pt = Rot(z) * Rot(y) * P,   rotate about y axes then z
    // Since we're undoing the rotation described in the definition
    // of the coordinate system,
    // we first rotate by ϑ = -(90 + spole_.lat()) around the y axis
    // (along the rotated Greenwich meridian)
    // and then by φ = -spole_.lon() degrees around the z axis):
    // (xt)   ( cos(φ), sin(φ), 0) (  cos(ϑ), 0, sin(ϑ)) (x)
    // (yt) = (-sin(φ), cos(φ), 0).(  0     , 1, 0     ).(y)
    // (zt)   ( 0     , 0     , 1) ( -sin(ϑ), 0, cos(ϑ)) (z)

    // Expanded
    // xt =  cos(ϑ) cos(φ) x + sin(φ) y + sin(ϑ) cos(φ) z
    // yt = -cos(ϑ) sin(φ) x + cos(φ) y - sin(ϑ) sin(φ) z
    // zt = -sin(ϑ)        x            + cos(ϑ)        z

    rotate_ = RotationMatrix{{{cos_theta * cos_phi, sin_phi, sin_theta * cos_phi},
                              {-cos_theta * sin_phi, cos_phi, -sin_theta * sin_phi},
                              {-sin_theta, 0., cos_theta}}};

    // Assume right hand rule, rotate about z axes and then y
    // P = Rot(y) * Rot(z) * Pt
    // x   (  cos(ϑ), 0, -sin(ϑ)) ( cos(φ), -sin(φ), 0) (xt)
    // y = (  0     , 1,  0     ) ( sin(φ), cos(φ),  0) (yt)
    // z   ( sin(ϑ), 0,   cos(ϑ)) ( 0     , 0     ,  1) (zt)

    // Expanded
    // x   ( cos(ϑ)cos(φ) , -cos(ϑ)sin(φ) , -sin(ϑ)) (xt)
    // y = ( sin(φ)       ,  cos(φ)       ,  0     ).(yt)
    // z   ( sin(ϑ) cos(φ), -sin(ϑ) sin(φ),  cos(ϑ)) (zt)

    unrotate_ = RotationMatrix{{{cos_theta * cos_phi, -cos_theta * sin_phi, -sin_theta},
                                {sin_phi, cos_phi, 0.},
                                {sin_theta * cos_phi, -sin_theta * sin_phi, cos_theta}}};
}

Rotation::Rotation(const PointLonLat& south_pole, double rotation_angle) {
    spole_ = south_pole;
    npole_ = PointLonLat(spole_.lon() - 180., -spole_.lat());
    if (npole_.lon() < 0.) {
        npole_.lon() += 360.;
    }
    angle_ = wrap_angle(rotation_angle);

    precompute();
}

Rotation::Rotation(const eckit::Parametrisation& p) {
    // get rotation angle
    p.get("rotation_angle", angle_);
    angle_ = wrap_angle(angle_);

    // get pole
    std::vector<double> pole(2);
    if (p.get("north_pole", pole)) {
        npole_ = PointLonLat(pole.data());
        spole_ = PointLonLat(npole_.lon() - 180., -npole_.lat());
        if (spole_.lon() < 0.) {
            spole_.lon() += 360.;
        }
    }
    else if (p.get("south_pole", pole)) {
        spole_ = PointLonLat(pole.data());
        npole_ = PointLonLat(spole_.lon() - 180., -spole_.lat());
        if (npole_.lon() < 0.) {
            npole_.lon() += 360.;
        }
    }

    precompute();
}

using RotationMatrix = std::array<std::array<double, 3>, 3>;

inline PointXYZ rotate_geocentric(const PointXYZ& p, const RotationMatrix& R) {
    return PointXYZ(R[XX][XX] * p.x() + R[XX][YY] * p.y() + R[XX][ZZ] * p.z(),
                    R[YY][XX] * p.x() + R[YY][YY] * p.y() + R[YY][ZZ] * p.z(),
                    R[ZZ][XX] * p.x() + R[ZZ][YY] * p.y() + R[ZZ][ZZ] * p.z());
}

void Rotation::rotate(double crd[]) const {
    if (!rotated_) {
        return;
    }

    crd[LON] -= angle_;

    if (!rotation_angle_only_) {
        const PointLonLat L(wrap_latitude({crd[LON], crd[LAT]}));
        PointXYZ P;
        UnitSphere::convertSphericalToCartesian(L, P);

        const PointXYZ Pt = rotate_geocentric(P, rotate_);
        PointLonLat Lt;
        UnitSphere::convertCartesianToSpherical(Pt, Lt);

        crd[LON] = Lt.lon();
        crd[LAT] = Lt.lat();
    }
}

void Rotation::unrotate(double crd[]) const {
    if (!rotated_) {
        return;
    }

    if (!rotation_angle_only_) {
        const PointLonLat Lt(crd);
        PointXYZ Pt;
        UnitSphere::convertSphericalToCartesian(Lt, Pt);

        const PointXYZ P = rotate_geocentric(Pt, unrotate_);
        PointLonLat L;
        UnitSphere::convertCartesianToSpherical(P, L);

        crd[LON] = L.lon();
        crd[LAT] = L.lat();
    }

    crd[LON] += angle_;
}

}  // namespace util
}  // namespace atlas
