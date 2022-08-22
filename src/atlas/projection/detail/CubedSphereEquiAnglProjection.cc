/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiAnglProjection.h"

#include <cmath>
#include <iostream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace {

static constexpr bool debug     = false;  // constexpr so compiler can optimize `if ( debug ) { ... }` out
static constexpr double deg2rad = atlas::util::Constants::degreesToRadians();
static constexpr double rad2deg = atlas::util::Constants::radiansToDegrees();

// Define small number relative to 360.
constexpr double epsilon = std::numeric_limits<double>::epsilon() * 360.;

// Define some "fuzzy" comparison operators.

// a is approximately equal to b.
bool equal(double a, double b) {
    return std::abs(a - b) <= epsilon;
}
// a is less than b.
bool lessThan(double a, double b) {
    return a < b && !equal(a, b);
}
// a is greater than b.
bool greaterThan(double a, double b) {
    return a > b && !equal(a, b);
}
// a is less than or approximately equal to b.
bool lessEqual(double a, double b) {
    return a < b || equal(a, b);
}
// a is greater than or approximately equal to b.
bool greaterEqual(double a, double b) {
    return a > b || equal(a, b);
}


}  // namespace

namespace atlas {
namespace projection {
namespace detail {


// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::CubedSphereEquiAnglProjection(const eckit::Parametrisation& params):
    CubedSphereProjectionBase(params) {}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::xy2alphabeta(double crd[], idx_t t) const {
    // Get tile centre.
    const auto& xyCentre = getCubedSphereTiles().tileCentre(static_cast<size_t>(t));

    // Check that xy coordinate is within valid "+" shaped halo region.
    const auto inCross = [&](const double crd[]) -> bool {
        return (greaterEqual(crd[XX], xyCentre[XX] - 45.) && lessEqual(crd[XX], xyCentre[XX] + 45.)) ||
               (greaterEqual(crd[YY], xyCentre[YY] - 45.) && lessEqual(crd[YY], xyCentre[YY] + 45.));
    };
    if (!inCross(crd)) {
        std::stringstream sStream;
        sStream << "xy coordinate (" << crd[0] << ", " << crd[1] << ") is not in range for tile " << t << ".";
        ATLAS_THROW_EXCEPTION(sStream.str());
    }

    // Get alphaBeta Jacobian.
    const auto alphabetaJacobian = getCubedSphereTiles().tileJacobian(static_cast<size_t>(t)).inverse();

    // Set (alpha, beta) coord.
    const Point2 alphabeta = alphabetaJacobian * (Point2(crd) - xyCentre);
    crd[0]                 = alphabeta[0];
    crd[1]                 = alphabeta[1];

    // Define correction.
    const auto correction = [](const double crd[]) -> double {
        return rad2deg * std::atan(std::tan(crd[0] * deg2rad) * std::tan(crd[1] * deg2rad));
    };

    // Correct halo (alpha, beta) coord.
    if (lessThan(crd[0], -45.)) {
        // Left.
        crd[1] = -correction(crd);
    }
    else if (greaterThan(crd[0], 45.)) {
        // Right.
        crd[1] = correction(crd);
    }
    else if (lessThan(crd[1], -45.)) {
        // Bottom.
        crd[0] = -correction(crd);
    }
    else if (greaterThan(crd[1], 45.)) {
        // Top.
        crd[0] = correction(crd);
    }
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::alphabeta2xy(double crd[], idx_t t) const {
    // Define correction.
    const auto correction1 = [](const double crd[]) -> double {
        return rad2deg * std::atan(std::tan(crd[1] * deg2rad) / std::tan(crd[0] * deg2rad));
    };
    const auto correction2 = [](const double crd[]) -> double {
        return rad2deg * std::atan(std::tan(crd[0] * deg2rad) / std::tan(crd[1] * deg2rad));
    };

    // Correct halo (alpha, beta) coord.
    if (lessThan(crd[0], -45.) && greaterThan(crd[1], crd[0]) && lessEqual(crd[1], -crd[0])) {
        // Left trapezium.
        crd[1] = -correction1(crd);
    }
    else if (greaterThan(crd[0], 45.) && greaterEqual(crd[1], -crd[0]) && lessThan(crd[1], crd[0])) {
        // Right trapezium.
        crd[1] = correction1(crd);
    }
    else if (lessThan(crd[1], -45.) && greaterEqual(crd[0], crd[1]) && lessThan(crd[0], -crd[1])) {
        // Bottom trapezium.
        crd[0] = -correction2(crd);
    }
    else if (greaterThan(crd[1], 45.) && greaterThan(crd[0], -crd[1]) && lessEqual(crd[0], crd[1])) {
        // Top trapezium.
        crd[0] = correction2(crd);
    }

    // Get tile centre.
    const auto& xyCentre = getCubedSphereTiles().tileCentre(static_cast<size_t>(t));

    // Get xy Jacobian.
    const auto xyJacobian = getCubedSphereTiles().tileJacobian(static_cast<size_t>(t));

    // Set xy coord.
    const Point2 xy = xyJacobian * Point2(crd) + xyCentre;
    crd[XX]         = xy[XX];
    crd[YY]         = xy[YY];
}

// -------------------------------------------------------------------------------------------------

Jacobian CubedSphereEquiAnglProjection::jacobian(const PointLonLat& lonlat, idx_t t) const {
    // Note: angular units cancel, so we leave all values in radians.

    // Convert lonlat to xyz on unit sphere.
    const double lambda = lonlat.lon() * deg2rad;
    const double phi    = lonlat.lat() * deg2rad;
    auto xyz            = PointXYZ{std::cos(lambda) * std::cos(phi), std::sin(lambda) * std::cos(phi), -std::sin(phi)};

    // Get derivatives of xyz with respect to lambda and phi.
    auto dxyz_by_dlambda = PointXYZ{-std::sin(lambda) * std::cos(phi), std::cos(lambda) * std::cos(phi), 0.};
    auto dxyz_by_dphi = PointXYZ{-std::cos(lambda) * std::sin(phi), -std::sin(lambda) * std::sin(phi), -std::cos(phi)};

    // Rotate vectors.
    const auto& tiles = getCubedSphereTiles();
    tiles.unrotate(t, xyz.data());
    tiles.unrotate(t, dxyz_by_dlambda.data());
    tiles.unrotate(t, dxyz_by_dphi.data());

    // Get derivatives of a and b with respect to xyz.
    // Note: a and b are xy displacements from the tile centre, *not* the
    // (alpha, beta) coordinates defined in the tileJacobian method.

    const double inv_x2y2 = 1. / (xyz.x() * xyz.x() + xyz.y() * xyz.y());
    const double inv_x2z2 = 1. / (xyz.x() * xyz.x() + xyz.z() * xyz.z());

    const auto da_by_dxyz = PointXYZ{-xyz.y() * inv_x2y2, xyz.x() * inv_x2y2, 0.};
    const auto db_by_dxyz = PointXYZ{xyz.z() * inv_x2z2, 0., -xyz.x() * inv_x2z2};

    // Use chain rule to get Jacobian.
    const auto& dot = eckit::geometry::Point3::dot;
    return Jacobian{{dot(da_by_dxyz, dxyz_by_dlambda), dot(da_by_dxyz, dxyz_by_dphi)},
                    {dot(db_by_dxyz, dxyz_by_dlambda), dot(db_by_dxyz, dxyz_by_dphi)}};
}

Jacobian CubedSphereEquiAnglProjection::alphabetaJacobian(const PointLonLat& lonlat, idx_t t) const {
    const auto& tiles = getCubedSphereTiles();
    return tiles.tileJacobian(t).inverse() * jacobian(lonlat, t);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::lonlat2xy(double crd[]) const {
    if (debug) {
        Log::info() << "equiangular lonlat2xy start : lonlat = " << crd[LON] << " " << crd[LAT] << std::endl;
    }

    idx_t t;
    double ab[2];   // alpha-beta coordinate
    double xyz[3];  // on Cartesian grid

    CubedSphereProjectionBase::lonlat2xy_pre(crd, t, xyz);

    // should be between -45.0 and 45.0
    // now calculate (alpha, beta) in radians.
    ab[0] = std::atan2(xyz[YY], xyz[XX]) * rad2deg;
    ab[1] = std::atan2(-xyz[ZZ], xyz[XX]) * rad2deg;  // I think the minus is here due to the
    // left coordinate system

    if (debug) {
        Log::info() << "equiangular lonlat2xy xyz ab : " << xyz[XX] << " " << xyz[YY] << " " << xyz[ZZ] << " "
                    << ab[LON] << " " << ab[LAT] << std::endl;
    }

    CubedSphereProjectionBase::alphabetat2xy(t, ab, crd);

    if (debug) {
        Log::info() << "equiangular lonlat2xy end : xy = " << crd[LON] << " " << crd[LAT] << std::endl;
    }
}

// -------------------------------------------------------------------------------------------------
// input should be Willems xy coordinate in degrees
//
void CubedSphereEquiAnglProjection::xy2lonlat(double crd[]) const {
    if (debug) {
        Log::info() << "xy2lonlat start xy = " << crd[LON] << " " << crd[LAT] << std::endl;
    }

    static const double rsq3 = 1.0 / std::sqrt(3.0);
    double xyz[3];
    double ab[2];  // alpha-beta coordinate
    idx_t t;       // tile index

    // calculate xy (in degrees) to alpha beta (in radians) and t - tile index.
    CubedSphereProjectionBase::xy2alphabetat(crd, t, ab);

    if (debug) {
        Log::debug() << "equiangular xy2lonlat:: crd t ab  : " << crd[0] << " " << crd[1] << " " << t << " " << ab[0]
                     << " " << ab[1] << std::endl;
    }

    xyz[0] = -rsq3;
    xyz[1] = -rsq3 * std::tan(ab[0] * deg2rad);
    xyz[2] = -rsq3 * std::tan(ab[1] * deg2rad);

    CubedSphereProjectionBase::xy2lonlat_post(xyz, t, crd);

    if (debug) {
        Log::info() << "end of equiangular xy2lonlat lonlat = " << crd[LON] << " " << crd[LAT] << std::endl;
    }
}

// -------------------------------------------------------------------------------------------------

Jacobian CubedSphereEquiAnglProjection::jacobian(const PointLonLat& lonlat) const {
    const auto& tiles = getCubedSphereTiles();
    const idx_t t     = tiles.indexFromLonLat(lonlat.data());
    return jacobian(lonlat, t);
}

// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::Spec CubedSphereEquiAnglProjection::spec() const {
    // Fill projection specification
    Spec proj;
    proj.set("type", static_type());
    return proj;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::hash(eckit::Hash& h) const {
    // Add to hash
    h.add(static_type());
    CubedSphereProjectionBase::hash(h);
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiAnglProjection> register_1(CubedSphereEquiAnglProjection::static_type());
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
