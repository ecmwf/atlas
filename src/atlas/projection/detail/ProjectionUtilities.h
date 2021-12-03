/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <functional>

#include "eckit/geometry/Sphere.h"

#include "atlas/domain.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

struct ProjectionUtilities {
    // -----------------------------------------------------------------------------------------------

    enum class CoordinateSystem
    {
        RIGHT_HAND,
        LEFT_HAND,
    };

    static constexpr bool debug = false;  // constexpr so compiler can optimize `if ( debug ) { ... }` out

    // -----------------------------------------------------------------------------------------------

    static void cartesianToSpherical(const double xyz[], double lonlat[], const CoordinateSystem coordinate_system,
                                     const double& radius = 0) {
        using eckit::geometry::Sphere;
        using util::Constants;

        // Make point objects.
        const auto pointXYZ = PointXYZ(xyz);
        auto pointLonLat    = PointLonLat();

        // Transform coordinates.
        auto r = radius != 0. ? radius : PointXYZ::norm(pointXYZ);
        Sphere::convertCartesianToSpherical(r, pointXYZ, pointLonLat);

        // Copy to array.
        lonlat[LON] = pointLonLat.lon();
        lonlat[LAT] = -pointLonLat.lat();

        // Left or right hand system.
        if (coordinate_system == CoordinateSystem::RIGHT_HAND) {
            lonlat[LAT] += 90.;
        }
    }

    //------------------------------------------------------------------------------------------------

    static void sphericalToCartesian(const double lonlat[], double xyz[], const CoordinateSystem coordinate_system,
                                     const double& radius = 0) {
        using eckit::geometry::Sphere;
        using util::Constants;

        // Make point objects.

        const auto pointLonLat = PointLonLat(lonlat);

        auto pointXYZ = PointXYZ();

        // Set Radius
        auto r = radius != 0 ? radius : util::Earth::radius();

        // Transform coordinates.
        Sphere::convertSphericalToCartesian(r, pointLonLat, pointXYZ, 0.0);

        if (debug) {
            Log::info() << "sphericalToCartesian:: pointLonLat pointXYZ = " << pointLonLat << " " << pointXYZ
                        << std::endl;
        }

        // Copy to array.
        xyz[XX] = pointXYZ.x();
        xyz[YY] = pointXYZ.y();
        xyz[ZZ] = pointXYZ.z();

        // Left or right hand system
        if (coordinate_system != CoordinateSystem::RIGHT_HAND) {
            xyz[ZZ] *= -1;
        }
    }

    //------------------------------------------------------------------------------------------------

    static void rotate3dX(const double angle, double xyz[]) {
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        double xyz_in[3];
        std::copy(xyz, xyz + 3, xyz_in);
        xyz[YY] = c * xyz_in[YY] + s * xyz_in[ZZ];
        xyz[ZZ] = -s * xyz_in[YY] + c * xyz_in[ZZ];
    };

    //------------------------------------------------------------------------------------------------

    static void rotate3dY(const double angle, double xyz[]) {
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        double xyz_in[3];
        std::copy(xyz, xyz + 3, xyz_in);
        xyz[XX] = c * xyz_in[XX] - s * xyz_in[ZZ];
        xyz[ZZ] = s * xyz_in[XX] + c * xyz_in[ZZ];
    };

    //------------------------------------------------------------------------------------------------

    static void rotate3dZ(const double angle, double xyz[]) {
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        double xyz_in[3];
        std::copy(xyz, xyz + 3, xyz_in);
        xyz[XX] = c * xyz_in[XX] + s * xyz_in[YY];
        xyz[YY] = -s * xyz_in[XX] + c * xyz_in[YY];
    };

    //------------------------------------------------------------------------------------------------
};

//--------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
