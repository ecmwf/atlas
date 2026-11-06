/*
 * (C) British Crown Copyright 2026 Met Office
 *
 */

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

/// @brief Partitioner that distributes grid points along a Hilbert space-filling curve.
///
/// All grid points are sorted by their Hilbert-curve index (computed over the
/// grid bounding box) and the resulting one-dimensional ordering is cut into
/// @c nb_partitions contiguous segments. Because the Hilbert curve preserves
/// spatial locality, each segment corresponds to a compact, contiguous
/// sub-domain, which keeps halo/communication volume low.
///
/// The partition is computed identically and independently on every MPI rank
/// (no communication), so it is fully deterministic and reproducible.
///
/// By default each grid point carries unit weight, so the curve is cut into
/// segments of (near-)equal point count. A user-supplied weight function can
/// instead balance an arbitrary, spatially-varying cost per point (for example
/// observation density, variable-resolution work, or masked/inactive points):
/// the curve is then cut into contiguous segments of (near-)equal total weight.
/// The unit-weight default exactly reproduces the equal-count behaviour, so the
/// weighting is a strict, opt-in generalisation.
///
/// Configuration options:
///   - "recursion" : <int> (default=30)  // Hilbert-curve recursion depth. Must be
///                                        // large enough to give each grid point a
///                                        // unique index (30 => 64-bit indices).
///
/// The weight function is not part of the (string/config-driven) factory
/// interface; it must be set programmatically via set_weight_function() (analytical
/// weights) or set_weights() (a precomputed per-point array) on a concrete instance.
/// This keeps the generic Partitioner API unchanged while still allowing callers
/// that know about weighting to opt in.
class SpaceFillingCurvePartitioner : public Partitioner {
public:
    /// @brief Weight function mapping a grid point (its global index in grid-iteration
    /// order, and its (lon,lat) coordinates in degrees) to a non-negative weight.
    /// Negative results are clamped to zero. The index makes it straightforward to
    /// look up a precomputed per-point array; the coordinate supports analytical
    /// weights. Either argument may be ignored.
    using WeightFunction = std::function<double(gidx_t index, const PointXY& lonlat)>;

    SpaceFillingCurvePartitioner();

    SpaceFillingCurvePartitioner(int N, const eckit::Parametrisation&);

    SpaceFillingCurvePartitioner(const eckit::Parametrisation&);

    using Partitioner::partition;
    void partition(const Grid&, int part[]) const override;

    std::string type() const override { return "hilbert"; }

    /// @brief Set a per-point weight function. When set, the curve is cut into
    /// segments of (near-)equal total weight rather than equal point count.
    /// Passing an empty function (the default) restores unit-weight behaviour.
    void set_weight_function(WeightFunction weight) { weight_ = std::move(weight); }

    /// @brief Set precomputed per-point weights, indexed in grid-iteration order
    /// (i.e. weights[n] is the weight of the n-th point of grid.lonlat()). The
    /// container must have exactly grid.size() entries; this is checked at
    /// partition() time. This is a convenience wrapper around set_weight_function()
    /// for the common case of a weight array (e.g. observation counts per point).
    void set_weights(std::vector<double> weights) {
        weight_ = [w = std::move(weights)](gidx_t index, const PointXY&) -> double { return w.at(index); };
    }

    /// @brief Whether a (non-empty) weight function has been set.
    bool has_weight_function() const { return static_cast<bool>(weight_); }

private:
    void setup(const eckit::Parametrisation&);

    idx_t recursion_{30};
    WeightFunction weight_;
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
