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
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace projection {
namespace detail {
class CubedSphereProjectionBase;
} // namespace detail
} // namespace projection
} // namespace atlas



namespace atlas {
namespace interpolation {
namespace method {
namespace cubedsphere {

using namespace util;

class CellFinder {
public:

    enum CellType {QUAD, TRIAG, INVALID};

    struct Cell {
        idx_t idx;
        Intersect isect;
        CellType type;
    };

    constexpr idx_t invalidIndex() const {return -1;}

    CellFinder(const Mesh& mesh, const util::Config& config = util::Config("include halo", false));

    Cell getCell(const PointXY& xy, double edgeEpsion = 5. * std::numeric_limits<double>::epsilon(),
                 double Epsion = 5. * std::numeric_limits<double>::epsilon()) const;

    Cell getCell(const PointLonLat& lonlat, double edgeEpsion = 5. * std::numeric_limits<double>::epsilon(),
                 double Epsion = 5. * std::numeric_limits<double>::epsilon()) const;




private:

    Mesh mesh_{};
    const projection::detail::CubedSphereProjectionBase* projection_{};
    std::array<util::IndexKDTree2D, 6> trees_{};



};

} // namespace cubedsphere
} // namespace method
} // namespace interpolation
} // namespace atlas
