/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereProjectionBase.h"

#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>

#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------
// Helper functions and variables local to this translation unit
namespace {

static constexpr double deg2rad = util::Constants::degreesToRadians();

static void schmidtTransform(double stretchFac, double targetLon, double targetLat, double lonlat[]) {
    double c2p1 = 1.0 + stretchFac * stretchFac;
    double c2m1 = 1.0 - stretchFac * stretchFac;

    double sin_p = std::sin(targetLat * deg2rad);
    double cos_p = std::cos(targetLat * deg2rad);

    double sin_lat;
    double cos_lat;
    double lat_t;

    if (std::abs(c2m1) > 1.0e-7) {
        sin_lat = std::sin(lonlat[LAT] * deg2rad);
        lat_t   = std::asin((c2m1 + c2p1 * sin_lat) / (c2p1 + c2m1 * sin_lat));
    }
    else {  // no stretching
        lat_t = lonlat[LAT];
    }

    sin_lat      = std::sin(lat_t);
    cos_lat      = std::cos(lat_t);
    double sin_o = -(sin_p * sin_lat + cos_p * cos_lat * cos(lonlat[LON] * deg2rad));

    if ((1. - std::abs(sin_o)) < 1.0e-7) {  // poles
        lonlat[LON] = 0.0;
        lonlat[LAT] = std::copysign(90.0, sin_o);
    }
    else {
        lonlat[LAT] = std::asin(sin_o);
        lonlat[LON] = targetLon + atan2(-cos_lat * std::sin(lonlat[LON] * deg2rad),
                                        -sin_lat * cos_p + cos_lat * sin_p * std::cos(lonlat[LON] * deg2rad));
        if (lonlat[LON] < 0.0) {
            lonlat[LON] += 360.;
        }
        else if (lonlat[LON] >= 360.) {
            lonlat[LON] -= 360.;
        }
    }
}

void sphericalToCartesian(const double lonlat[], double xyz[]) {
    auto crd_sys            = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 1.;
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz, crd_sys, radius);
}

void cartesianToSpherical(const double xyz[], double lonlat[]) {
    auto crd_sys            = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 0.;  // --> equivalent to radius = norm(xyz)
    ProjectionUtilities::cartesianToSpherical(xyz, lonlat, crd_sys, radius);
}

std::string getTileType(const eckit::Parametrisation& params) {
    std::string tileStr;
    params.get("tile.type", tileStr);
    return tileStr;
}

}  // namespace


// -------------------------------------------------------------------------------------------------

CubedSphereProjectionBase::CubedSphereProjectionBase(const eckit::Parametrisation& params):
    tiles_(getTileType(params)),
    tiles_offsets_ab2xy_(tiles_.ab2xyOffsets()),
    tiles_offsets_xy2ab_(tiles_.xy2abOffsets()) {
    ATLAS_TRACE("CubedSphereProjectionBase::CubedSphereProjectionBase");

    // Shift projection by a longitude
    shiftLon_ = 0.0;
    if (params.has("ShiftLon")) {
        params.get("ShiftLon", shiftLon_);
        ATLAS_ASSERT(shiftLon_ <= 90.0, "ShiftLon should be <= 90.0 degrees");
        ATLAS_ASSERT(shiftLon_ >= -90.0, "ShiftLon should be >= -90.0 degrees");
    }

    // Apply a Schmidt transform
    doSchmidt_  = false;
    stretchFac_ = 0.0;
    targetLon_  = 0.0;
    targetLat_  = 0.0;
    if (params.has("DoSchmidt")) {
        params.get("DoSchmidt", doSchmidt_);
        if (doSchmidt_) {
            params.get("StretchFac", stretchFac_);
            params.get("TargetLon", targetLon_);
            params.get("TargetLat", targetLat_);
        }
    }
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::hash(eckit::Hash& h) const {
    // Add stretching options to hash
    h.add(shiftLon_);
    h.add(doSchmidt_);
    h.add(stretchFac_);
    h.add(targetLon_);
    h.add(targetLat_);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat_post(double xyz[], const idx_t t, double crd[]) const {
    cartesianToSpherical(xyz, crd);

    if (crd[LON] < 0.) {
        crd[LON] += 360.;
    }
    crd[LON] -= 180.;

    // Convert to cartesian
    sphericalToCartesian(crd, xyz);

    // Perform tile specific rotation
    tiles_.rotate(t, xyz);

    // Back to lonlat
    cartesianToSpherical(xyz, crd);

    // Shift longitude
    if (shiftLon_ != 0.0) {
        crd[LON] = crd[LON] + shiftLon_;
        if (crd[LON] < -180.) {
            crd[LON] = 360. + crd[LON];
        }
        if (crd[LON] > 180.) {
            crd[LON] = -360. + crd[LON];
        }
    }

    // To 0, 360
    if (crd[LON] < 0.0) {
        crd[LON] += 360.;
    }

    // Schmidt transform
    if (doSchmidt_) {
        schmidtTransform(stretchFac_, targetLon_, targetLat_, crd);
    }

    // longitude does not make sense at the poles - set to 0.
    if (std::abs(std::abs(crd[LAT]) - 90.) < 1e-15) {
        crd[LON] = 0.;
    }
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xy_pre(double crd[], idx_t& t, double xyz[]) const {
    if (std::abs(crd[LON]) < 1e-15) {
        crd[LON] = 0.;
    }
    if (std::abs(crd[LAT]) < 1e-15) {
        crd[LAT] = 0.;
    }

    // To [-45.0, 315)
    if (crd[LON] >= 315.0) {
        crd[LON] -= 360.;
    }

    // find tile which this lonlat is linked to
    // works [-45, 315.0)
    t = tiles_.indexFromLonLat(crd);

    sphericalToCartesian(crd, xyz);
    tiles_.unrotate(t, xyz);
}

// -------------------------------------------------------------------------------------------------


void CubedSphereProjectionBase::xy2alphabetat(const double xy[], idx_t& t, double ab[]) const {
    // xy is in degrees while ab is in degree
    // ab are the  (alpha, beta) coordinates and t is the tile index.

    t                  = tiles_.indexFromXY(xy);
    double normalisedX = xy[XX] / 90.;
    double normalisedY = (xy[YY] + 135.) / 90.;
    ab[LON]            = (normalisedX - tiles_offsets_xy2ab_[XX][size_t(t)]) * 90.0 - 45.0;
    ab[LAT]            = (normalisedY - tiles_offsets_xy2ab_[YY][size_t(t)]) * 90.0 - 45.0;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::alphabetat2xy(const idx_t t, const double ab[], double xy[]) const {
    // xy and ab are in degrees
    // (alpha, beta) and tiles.

    xy[XX] = ab[LON] + 45.0 + tiles_offsets_ab2xy_[LON][size_t(t)];
    xy[YY] = ab[LAT] + 45.0 + tiles_offsets_ab2xy_[LAT][size_t(t)];

    tiles_.enforceXYdomain(xy);
}

// -------------------------------------------------------------------------------------------------

Point2 CubedSphereProjectionBase::xy2alphabeta(const Point2& xy, idx_t t) const {
    auto alphabeta = Point2(xy);
    xy2alphabeta(alphabeta.data(), t);
    return alphabeta;
}

// -------------------------------------------------------------------------------------------------

Point2 CubedSphereProjectionBase::alphabeta2xy(const Point2& alphabeta, idx_t t) const {
    auto xy = Point2(alphabeta);
    alphabeta2xy(xy.data(), t);
    return xy;
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
