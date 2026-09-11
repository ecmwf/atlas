/*
 * (C) British Crown Copyright 2026 Met Office
 *
 */

#include "atlas/grid/detail/partitioner/SpaceFillingCurvePartitioner.h"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "atlas/domain/Domain.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/SpaceFillingCurve.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

SpaceFillingCurvePartitioner::SpaceFillingCurvePartitioner(): Partitioner() {}

SpaceFillingCurvePartitioner::SpaceFillingCurvePartitioner(int N, const eckit::Parametrisation& config):
    Partitioner(N, config) {
    setup(config);
}

SpaceFillingCurvePartitioner::SpaceFillingCurvePartitioner(const eckit::Parametrisation& config): Partitioner(config) {
    setup(config);
}

void SpaceFillingCurvePartitioner::setup(const eckit::Parametrisation& config) {
    config.get("recursion", recursion_);
    ATLAS_ASSERT(recursion_ > 0, "SpaceFillingCurvePartitioner: 'recursion' must be strictly positive");
}

void SpaceFillingCurvePartitioner::partition(const Grid& grid, int part[]) const {
    const int nb_parts     = static_cast<int>(nb_partitions());
    const gidx_t nb_points = grid.size();

    if (nb_parts == 1) {  // trivial solution, so much faster
        for (gidx_t j = 0; j < nb_points; ++j) {
            part[j] = 0;
        }
        return;
    }

    ATLAS_TRACE("SpaceFillingCurvePartitioner::partition");

    ATLAS_ASSERT(nb_parts <= nb_points,
                 "SpaceFillingCurvePartitioner: cannot partition a grid into more parts than it has points");

    // Gather all point coordinates and compute the global bounding box in a single
    // pass. This is done locally and identically on every MPI rank, so the resulting
    // partition is deterministic and reproducible regardless of rank count.
    std::vector<PointXY> points;
    points.reserve(nb_points);
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
    for (const auto& p : grid.lonlat()) {
        points.emplace_back(p[0], p[1]);
        xmin = std::min(xmin, p[0]);
        xmax = std::max(xmax, p[0]);
        ymin = std::min(ymin, p[1]);
        ymax = std::max(ymax, p[1]);
    }
    // Guard against a degenerate (zero-width) bounding box.
    if (!(xmax > xmin)) {
        xmax = xmin + 1.;
    }
    if (!(ymax > ymin)) {
        ymax = ymin + 1.;
    }

    HilbertCurve hilbert{RectangularDomain{{xmin, xmax}, {ymin, ymax}}, recursion_};

    // (hilbert key, original grid index)
    std::vector<std::pair<gidx_t, gidx_t>> sfc;
    sfc.reserve(nb_points);
    ATLAS_TRACE_SCOPE("compute hilbert keys") {
        for (gidx_t n = 0; n < nb_points; ++n) {
            sfc.emplace_back(hilbert(points[n]), n);
        }
    }

    ATLAS_TRACE_SCOPE("sort") {
        std::sort(sfc.begin(), sfc.end());
    }

    // Per-point weights (default: unit weight => equal-count cut). When a weight
    // function is supplied, negative weights are clamped to zero; if the total
    // weight is non-positive (e.g. every point masked out) we fall back to unit
    // weights so the partition is still well-defined.
    std::vector<double> weights;
    double total_weight = static_cast<double>(nb_points);
    if (weight_) {
        weights.resize(nb_points);
        total_weight = 0.;
        for (gidx_t n = 0; n < nb_points; ++n) {
            double w = weight_(n, points[n]);
            if (!(w > 0.)) {
                w = 0.;
            }
            weights[n] = w;
            total_weight += w;
        }
        if (!(total_weight > 0.)) {
            std::fill(weights.begin(), weights.end(), 1.);
            total_weight = static_cast<double>(nb_points);
        }
    }

    // Cut the sorted curve into nb_parts contiguous segments of (near-)equal total
    // weight. The boundary test is written in the scaled form
    //     accumulated * nb_parts >= (p+1) * total_weight
    // rather than accumulating a running target, so that it introduces no
    // floating-point drift. In particular, with unit weights the accumulated value
    // and total_weight are exact integers, making this reduce *exactly* to an
    // equal-count partition where each partition receives floor(N/P) or ceil(N/P)
    // points. Each point is assigned to the current partition first; once the
    // accumulated weight reaches the target, subsequent points advance to the next
    // partition, so the overshoot at every boundary is bounded by one point's weight.
    double accumulated = 0.;
    int p              = 0;
    for (gidx_t k = 0; k < nb_points; ++k) {
        const gidx_t index = sfc[k].second;
        part[index]        = p;
        accumulated += weights.empty() ? 1. : weights[index];
        while (p < nb_parts - 1 && accumulated * nb_parts >= static_cast<double>(p + 1) * total_weight) {
            ++p;
        }
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::SpaceFillingCurvePartitioner>
    __SpaceFillingCurve("hilbert");
}
