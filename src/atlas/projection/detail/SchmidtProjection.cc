/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#include <cmath>

#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/SchmidtProjection.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/UnitSphere.h"

namespace {
static double D2R(const double x) {
    return atlas::util::Constants::degreesToRadians() * x;
}
static double R2D(const double x) {
    return atlas::util::Constants::radiansToDegrees() * x;
}
}  // namespace

namespace atlas {
namespace projection {
namespace detail {

PointLonLat rotation_north_pole(const Rotated& rotation) {
    return rotation.northPole();
}

PointLonLat rotation_north_pole(const NotRotated& rotation) {
    return PointLonLat{0., 90.};
}

// constructor
template <typename Rotation>
SchmidtProjectionT<Rotation>::SchmidtProjectionT(const eckit::Parametrisation& params):
    ProjectionImpl(), rotation_(params) {
    if (!params.get("stretching_factor", c_)) {
        throw_Exception("stretching_factor missing in Params", Here());
    }
    ATLAS_ASSERT(c_ != 0.);

    north0_ = {0.0, 0.0, 1.0};

    atlas::util::UnitSphere::convertSphericalToCartesian(rotation_north_pole(rotation_), north1_);
    north1_ = PointXYZ::normalize(north1_);
}

// constructor
template <typename Rotation>
SchmidtProjectionT<Rotation>::SchmidtProjectionT(): ProjectionImpl(), rotation_(util::NoConfig()) {}

template <typename Rotation>
void SchmidtProjectionT<Rotation>::xy2lonlat(double crd[]) const {
    // stretch
    crd[1] = R2D(std::asin(std::cos(2. * std::atan(1 / c_ * std::tan(std::acos(std::sin(D2R(crd[1]))) * 0.5)))));

    // perform rotation
    rotation_.rotate(crd);
}

template <typename Rotation>
void SchmidtProjectionT<Rotation>::lonlat2xy(double crd[]) const {
    // inverse rotation
    rotation_.unrotate(crd);

    // unstretch
    crd[1] = R2D(std::asin(std::cos(2. * std::atan(c_ * std::tan(std::acos(std::sin(D2R(crd[1]))) * 0.5)))));
}

template <>
ProjectionImpl::Jacobian SchmidtProjectionT<NotRotated>::jacobian(const PointLonLat&) const {
    throw_NotImplemented("SchmidtProjectionT<NotRotated>::jacobian", Here());
}

template <typename Rotation>
ProjectionImpl::Jacobian SchmidtProjectionT<Rotation>::jacobian(const PointLonLat& lonlat) const {
    double xy[2] = {lonlat.lon(), lonlat.lat()};

    lonlat2xy(xy);

    PointXYZ xyz;
    atlas::util::UnitSphere::convertSphericalToCartesian(lonlat, xyz);


    double zomc2 = 1.0 - 1.0 / (c_ * c_);
    double zopc2 = 1.0 + 1.0 / (c_ * c_);

    double zcosy   = std::cos(D2R(xy[1]));
    double zsiny   = std::sin(D2R(xy[1]));
    double zcoslat = std::cos(D2R(lonlat.lat()));

    double zfactor = std::sqrt((zopc2 + zsiny * zomc2) * (zopc2 + zsiny * zomc2) / (zopc2 * zopc2 - zomc2 * zomc2));

    // Base vectors in unrotated frame
    auto u0 = PointXYZ::normalize(PointXYZ::cross(north0_, xyz));
    auto v0 = PointXYZ::normalize(PointXYZ::cross(xyz, u0));

    // Base vectors in rotated frame
    auto u1 = PointXYZ::normalize(PointXYZ::cross(north1_, xyz));
    auto v1 = PointXYZ::normalize(PointXYZ::cross(xyz, u1));

    double u0u1 = PointXYZ::dot(u0, u1);
    double v0u1 = PointXYZ::dot(v0, u1);

    double u0v1 = PointXYZ::dot(u0, v1);
    double v0v1 = PointXYZ::dot(v0, v1);

    Jacobian jac;
    jac[0] = {zcoslat * u0u1 * zfactor / zcosy, v0u1 * zfactor / zcosy};
    jac[1] = {zcoslat * u0v1 * zfactor, v0v1 * zfactor};

    return jac;
}

// specification
template <typename Rotation>
typename SchmidtProjectionT<Rotation>::Spec SchmidtProjectionT<Rotation>::spec() const {
    Spec proj_spec;
    proj_spec.set("type", static_type());
    proj_spec.set("stretching_factor", c_);
    rotation_.spec(proj_spec);
    return proj_spec;
}

template <typename Rotation>
void SchmidtProjectionT<Rotation>::hash(eckit::Hash& hsh) const {
    hsh.add(static_type());
    rotation_.hash(hsh);
    hsh.add(c_);
}

template class SchmidtProjectionT<NotRotated>;
template class SchmidtProjectionT<Rotated>;

namespace {
static ProjectionBuilder<SchmidtProjection> register_1(SchmidtProjection::static_type());
static ProjectionBuilder<RotatedSchmidtProjection> register_2(RotatedSchmidtProjection::static_type());
}  // namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
