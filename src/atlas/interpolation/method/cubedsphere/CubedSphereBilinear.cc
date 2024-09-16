/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/cubedsphere/CubedSphereBilinear.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/cubedsphere/CellFinder.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {
MethodBuilder<CubedSphereBilinear> __builder("cubedsphere-bilinear");
}  // namespace

void CubedSphereBilinear::do_setup(const Grid& source, const Grid& target, const Cache&) {
    ATLAS_NOTIMPLEMENTED;
}

void CubedSphereBilinear::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    source_ = source;
    target_ = target;

    const auto ncSource = functionspace::NodeColumns(source);

    ATLAS_ASSERT(ncSource);
    ATLAS_ASSERT(target_);

    // Enable or disable halo exchange.
    this->allow_halo_exchange_ = halo_exchange_;


    // return early if no output points on this partition reserve is called on
    // the triplets but also during the sparseMatrix constructor. This won't
    // work for empty matrices
    if (target_.size() == 0) {
        return;
    }

    const auto finder = cubedsphere::CellFinder(ncSource.mesh(), util::Config("halo", halo_));

    // Numeric tolerance should scale with N.
    const auto N         = CubedSphereGrid(ncSource.mesh().grid()).N();
    const auto tolerance = 2. * std::numeric_limits<double>::epsilon() * N;

    // Loop over target at calculate interpolation weights.
    auto weights          = std::vector<Triplet>{};
    const auto ghostView  = array::make_view<int, 1>(target_.ghost());
    const auto lonlatView = array::make_view<double, 2>(target_.lonlat());
    const auto tijView    = array::make_view<idx_t, 2>(ncSource.mesh().cells().field("tij"));

    // Make vector of tile indices for each target point (needed for vector field interpolation).
    std::vector<idx_t> tileIndex{};

    for (idx_t i = 0; i < target_.size(); ++i) {
        if (!ghostView(i)) {
            const auto cell =
                finder.getCell(PointLonLat(lonlatView(i, LON), lonlatView(i, LAT)), listSize_, tolerance, tolerance);

            if (!cell.isect) {
                ATLAS_THROW_EXCEPTION(
                    "Cannot find a cell surrounding target"
                    "point " +
                    std::to_string(i) + ".");
            }

            tileIndex.push_back(tijView(cell.idx, 0));
            const auto& isect = cell.isect;
            const auto& j     = cell.nodes;

            switch (cell.nodes.size()) {
                case (3): {
                    // Cell is a triangle.
                    weights.emplace_back(i, j[0], 1. - isect.u - isect.v);
                    weights.emplace_back(i, j[1], isect.u);
                    weights.emplace_back(i, j[2], isect.v);
                    break;
                }
                case (4): {
                    // Cell is quad.
                    weights.emplace_back(i, j[0], (1. - isect.u) * (1. - isect.v));
                    weights.emplace_back(i, j[1], isect.u * (1. - isect.v));
                    weights.emplace_back(i, j[2], isect.u * isect.v);
                    weights.emplace_back(i, j[3], (1. - isect.u) * isect.v);
                    break;
                }
                default: {
                    ATLAS_THROW_EXCEPTION("Unknown cell type with " + std::to_string(cell.nodes.size()) + " nodes.");
                }
            }
        }
    }

    // fill sparse matrix and return.
    Matrix A(target_.size(), source_.size(), weights);
    setMatrix(A);

    // Add tile index metadata to target.
    target_->metadata().set("tile index", tileIndex);
}

void CubedSphereBilinear::print(std::ostream&) const {
    ATLAS_NOTIMPLEMENTED;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
