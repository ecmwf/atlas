/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/knn/GridBoxMethod.h"

#include <algorithm>
#include <vector>

#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"


namespace atlas {
namespace interpolation {
namespace method {


GridBoxMethod::GridBoxMethod(const Method::Config& config): KNearestNeighboursBase(config) {
    config.get("matrix_free", matrixFree_ = false);
    config.get("fail_early", failEarly_ = true);
    config.get("gaussian_weighted_latitudes", gaussianWeightedLatitudes_ = true);
}


GridBoxMethod::~GridBoxMethod() = default;


void GridBoxMethod::print(std::ostream& out) const {
    out << "GridBoxMethod[]";
}

namespace {
StructuredGrid extract_structured_grid(const FunctionSpace& fs) {
    if (functionspace::StructuredColumns{fs}) {
        return functionspace::StructuredColumns(fs).grid();
    }
    if (functionspace::NodeColumns{fs}) {
        return functionspace::NodeColumns(fs).mesh().grid();
    }
    return Grid();
}
}

void GridBoxMethod::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
    if( mpi::size() > 1 ) {
        ATLAS_THROW_EXCEPTION("Cannot use GridBoxMethod in parallel yet.");
    }
    auto src_grid = extract_structured_grid(source);
    auto tgt_grid = extract_structured_grid(target);
    if( not src_grid ) {
        ATLAS_THROW_EXCEPTION("Could not extract StructuredGrid from source function space " << source.type() );
    }
    if( not tgt_grid ) {
        ATLAS_THROW_EXCEPTION("Could not extract StructuredGrid from target function space " << target.type() );
    }
    do_setup(src_grid,tgt_grid,Cache());
}

bool GridBoxMethod::intersect(size_t i, const GridBox& box, const util::IndexKDTree::ValueList& closest,
                              std::vector<eckit::linalg::Triplet>& triplets) const {
    ASSERT(!closest.empty());

    triplets.clear();
    triplets.reserve(closest.size());

    double area = box.area();
    ASSERT(area > 0.);

    double sumSmallAreas = 0.;
    for (auto& c : closest) {
        auto j = c.payload();
        ASSERT(j >= 0);

        auto smallBox = sourceBoxes_.at(size_t(j));
        if (box.intersects(smallBox)) {
            double smallArea = smallBox.area();
            ASSERT(smallArea > 0.);

            triplets.emplace_back(i, j, smallArea / area);
            sumSmallAreas += smallArea;

            if (eckit::types::is_approximately_equal(area, sumSmallAreas, 1. /*m^2*/)) {
                return true;
            }
        }
    }

    if (failEarly_) {
        Log::error() << "Failed to intersect grid box " << i << ", " << box << std::endl;
        throw_Exception("Failed to intersect grid box");
    }

    failures_.push_front(i);
    triplets.clear();
    return false;
}


void GridBoxMethod::do_setup(const Grid& source, const Grid& target, const Cache& cache) {
    ATLAS_TRACE("GridBoxMethod::setup()");

    if (mpi::size() > 1) {
        ATLAS_THROW_EXCEPTION("Not implemented for MPI-parallel runs");
    }

    ATLAS_ASSERT(source);
    ATLAS_ASSERT(target);

    functionspace::PointCloud src(source);
    functionspace::PointCloud tgt(target);
    ATLAS_ASSERT(src);
    ATLAS_ASSERT(tgt);
    source_ = src;
    target_ = tgt;

    if (not matrixFree_ && interpolation::MatrixCache(cache)) {
        setMatrix(cache);
        ATLAS_ASSERT(matrix().rows() == target.size());
        ATLAS_ASSERT(matrix().cols() == source.size());
        return;
    }

    if (not extractTreeFromCache(cache)) {
        buildPointSearchTree(src);
    }

    sourceBoxes_ = GridBoxes(source, gaussianWeightedLatitudes_);
    targetBoxes_ = GridBoxes(target, gaussianWeightedLatitudes_);

    searchRadius_ = sourceBoxes_.getLongestGridBoxDiagonal() + targetBoxes_.getLongestGridBoxDiagonal();
    failures_.clear();

    if (matrixFree_) {
        return;
    }

    std::vector<Triplet> allTriplets;

    {
        ATLAS_TRACE("GridBoxMethod::setup: intersecting grid boxes");

        constexpr double TIMED = 5.;
        eckit::ProgressTimer progress("Intersecting", targetBoxes_.size(), "grid box", TIMED);

        std::vector<Triplet> triplets;
        size_t i = 0;
        for (auto p : tgt.iterate().lonlat()) {
            ++progress;

            if (intersect(i, targetBoxes_.at(i), pTree_.closestPointsWithinRadius(p, searchRadius_), triplets)) {
                std::copy(triplets.begin(), triplets.end(), std::back_inserter(allTriplets));
            }

            ++i;
        }

        if (!failures_.empty()) {
            giveUp(failures_);
        }
    }

    {
        ATLAS_TRACE("GridBoxMethod::setup: build interpolant matrix");
        setMatrix(targetBoxes_.size(), sourceBoxes_.size(), allTriplets);
    }
}


void GridBoxMethod::giveUp(const std::forward_list<size_t>& failures) {
    Log::warning() << "Failed to intersect grid boxes: ";

    constexpr int COUNTED = 10;

    int count = 0;
    auto sep  = "";
    for (const auto& f : failures) {
        if (count++ < COUNTED) {
            Log::warning() << sep << f;
            sep = ", ";
        }
    }
    Log::warning() << "... (" << eckit::Plural(count, "total failure") << std::endl;

    throw_Exception("Failed to intersect grid boxes");
}

Cache GridBoxMethod::createCache() const {
    Cache cache;
    cache.add(interpolation::IndexKDTreeCache(pTree_));
    if (not matrix().empty()) {
        cache.add(Method::createCache());
    }
    return cache;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
