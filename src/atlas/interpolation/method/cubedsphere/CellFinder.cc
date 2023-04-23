/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/cubedsphere/CellFinder.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/interpolation/Vector2D.h"
#include "atlas/interpolation/element/Quad2D.h"
#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/Topology.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace cubedsphere {

CellFinder::CellFinder(const Mesh& mesh, const util::Config& config): mesh_{mesh} {
    // Check mesh and get projection.
    const auto csGrid = CubedSphereGrid(mesh_.grid());
    ATLAS_ASSERT_MSG(csGrid, "cubedsphere::CellFinder requires a cubed sphere mesh.");
    projection_ = &(csGrid.cubedSphereProjection());

    // Get views to cell data.
    const auto lonlatView   = array::make_view<double, 2>(mesh_.cells().field("lonlat"));
    const auto haloView = array::make_view<int, 1>(mesh_.cells().halo());

    // make points and payloads vectors.
    auto points   = std::vector<PointLonLat>{};
    auto payloads = std::vector<idx_t>{};

    // Iterate over cells.
    auto halo = config.getInt("halo", 0);
    for (idx_t i = 0; i < mesh_.cells().size(); ++i) {
        if (haloView(i) <= halo) {
            points.emplace_back(PointLonLat(lonlatView(i, LON), lonlatView(i, LAT)));
            payloads.emplace_back(i);
        }
    }

    // build cell kd-tree.
    tree_.build(points, payloads);
}

CellFinder::Cell CellFinder::getCell(const PointLonLat& lonlat, size_t listSize, double edgeEpsilon, double epsilon) const {
    // Convert xy to alphabeta and t;
    const auto& tiles    = projection_->getCubedSphereTiles();

    // Get mesh nodes and connectivity table.
    const auto nodeXyView            = array::make_view<double, 2>(mesh_.nodes().xy());
    const auto& nodeConnectivity     = mesh_.cells().node_connectivity();

    // Get four nearest cell-centres to xy.
    if (tree_.size() == 0) {
        return Cell{-1, {}, Intersect{}.fail()};
    }

    const auto valueList = tree_.closestPoints(lonlat, std::min(listSize, tree_.size()));

    // Get view to cell flags.
    const auto cellFlagsView = array::make_view<int, 1>(mesh_.cells().flags());
    const auto cellTijView   = array::make_view<idx_t, 2>(mesh_.cells().field("tij"));

    for (const auto& value : valueList) {
        const auto i = value.payload();

        // Ignore invalid cells.
        if ((cellFlagsView(i) & Topology::INVALID)) {
            break;
        }

        const auto t   = cellTijView(i,0);
        const auto row = nodeConnectivity.row(i);

        if (row.size() == 4) {
            auto quadAlphabeta = std::array<Vector2D, 4>{};
            for (size_t k = 0; k < 4; ++k) {
                const auto xyNode = PointXY(nodeXyView(row(k), XX), nodeXyView(row(k), YY));
                quadAlphabeta[k]  = Vector2D(projection_->xy2alphabeta(xyNode, t).data());
            }

            const auto quad =
                element::Quad2D(quadAlphabeta[0], quadAlphabeta[1], quadAlphabeta[2], quadAlphabeta[3]);

            const auto alphabeta = projection_->lonlat2alphabeta(lonlat, t);
            const auto isect = quad.localRemap(alphabeta, edgeEpsilon, epsilon);

            if (isect) {
                return Cell{i, {row(0), row(1), row(2), row(3)}, isect};
            }
        }
        else {
            // Cell is triangle.
            auto triagAlphabeta = std::array<Vector2D, 3>{};
            for (size_t k = 0; k < 3; ++k) {
                const auto xyNode = PointXY(nodeXyView(row(k), XX), nodeXyView(row(k), YY));
                triagAlphabeta[k] = Vector2D(projection_->xy2alphabeta(xyNode, t).data());
            }

            const auto triag = element::Triag2D(triagAlphabeta[0], triagAlphabeta[1], triagAlphabeta[2]);

            const auto alphabeta = projection_->lonlat2alphabeta(lonlat, t);
            const auto isect = triag.intersects(alphabeta, edgeEpsilon, epsilon);

            if (isect) {
                return Cell{i, {row(0), row(1), row(2)}, isect};
            }
        }
    }

    // Couldn't find a cell.
    return Cell{-1, {}, Intersect{}.fail()};
}

CellFinder::Cell CellFinder::getCell(const PointXY& xy, size_t listSize, double edgeEpsilon,
                                     double epsilon) const {
    const auto lonlat = projection_->lonlat(xy);
    return getCell(lonlat, listSize, edgeEpsilon, epsilon);
}


}  // namespace cubedsphere
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
