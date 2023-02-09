/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>

#include "atlas/interpolation/method/Intersect.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/Config.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace projection {
namespace detail {
class CubedSphereProjectionBase;
}  // namespace detail
}  // namespace projection
}  // namespace atlas


namespace atlas {
namespace interpolation {
namespace method {
namespace cubedsphere {

using namespace util;

/// @brief class to find points within cells of cubedsphere mesh.
class CellFinder {
public:
    struct Cell {
        idx_t idx;
        std::vector<idx_t> nodes;
        Intersect isect;
    };

    /// @brief Constructor.
    CellFinder(const Mesh& mesh, const util::Config& config = util::Config("halo", 0));

    /// @brief Find a cell which encompasses an xy point.
    Cell getCell(const PointXY& xy, size_t listSize = 4,
                 double edgeEpsilon = 5. * std::numeric_limits<double>::epsilon(),
                 double epsilon     = 5. * std::numeric_limits<double>::epsilon()) const;

    /// @brief Find a cell which encompasses a lonlat point.
    Cell getCell(const PointLonLat& lonlat, size_t listSize = 4,
                 double edgeEpsilon = 5. * std::numeric_limits<double>::epsilon(),
                 double epsilon     = 5. * std::numeric_limits<double>::epsilon()) const;


private:
    Mesh mesh_{};
    const projection::detail::CubedSphereProjectionBase* projection_{};
    util::IndexKDTree tree_{Geometry{}};
};

}  // namespace cubedsphere
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
