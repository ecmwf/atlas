/*
 * (C) British Crown Copyright 2026 Met Office
 *
 */

#include <algorithm>
#include <vector>

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/SpaceFillingCurvePartitioner.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using grid::detail::partitioner::SpaceFillingCurvePartitioner;

//----------------------------------------------------------------------------------------------------------------------

// Helper: every point is assigned to a valid partition in [0, nb_parts).
void check_valid_partitions(const grid::Distribution& d, int nb_parts) {
    for (gidx_t n = 0; n < d.size(); ++n) {
        const int p = d.partition(n);
        EXPECT(p >= 0);
        EXPECT(p < nb_parts);
    }
}

// Helper: partition sizes differ by at most 1 (equal-count balancing).
void check_balanced(const grid::Distribution& d, int nb_parts) {
    const auto& counts = d.nb_pts();
    EXPECT_EQ(static_cast<int>(counts.size()), nb_parts);

    const gidx_t size      = d.size();
    const idx_t floor_size = static_cast<idx_t>(size / nb_parts);
    const idx_t ceil_size  = floor_size + (size % nb_parts != 0 ? 1 : 0);

    idx_t total = 0;
    for (idx_t c : counts) {
        EXPECT(c >= floor_size);
        EXPECT(c <= ceil_size);
        total += c;
    }
    EXPECT_EQ(total, static_cast<idx_t>(size));
}

//----------------------------------------------------------------------------------------------------------------------

CASE("hilbert partitioner is registered") {
    EXPECT(grid::Partitioner::exists("hilbert"));
}

CASE("single partition assigns everything to rank 0") {
    Grid grid("O32");
    grid::Partitioner partitioner("hilbert", 1);
    auto d = partitioner.partition(grid);

    EXPECT_EQ(d.nb_partitions(), 1);
    for (gidx_t n = 0; n < d.size(); ++n) {
        EXPECT_EQ(d.partition(n), 0);
    }
}

CASE("all points assigned and partitions balanced") {
    Grid grid("O32");
    for (int nb_parts : {2, 3, 4, 7, 16}) {
        grid::Partitioner partitioner("hilbert", nb_parts);
        auto d = partitioner.partition(grid);

        EXPECT_EQ(d.nb_partitions(), nb_parts);
        EXPECT_EQ(d.size(), static_cast<gidx_t>(grid.size()));
        check_valid_partitions(d, nb_parts);
        check_balanced(d, nb_parts);
    }
}

CASE("partition is deterministic / reproducible") {
    Grid grid("O32");
    grid::Partitioner p1("hilbert", 5);
    grid::Partitioner p2("hilbert", 5);
    auto d1 = p1.partition(grid);
    auto d2 = p2.partition(grid);

    EXPECT_EQ(d1.size(), d2.size());
    for (gidx_t n = 0; n < d1.size(); ++n) {
        EXPECT_EQ(d1.partition(n), d2.partition(n));
    }
}

CASE("recursion configuration option is honoured") {
    Grid grid("O32");
    util::Config config;
    config.set("partitions", 4);
    config.set("recursion", 20);
    grid::Partitioner partitioner("hilbert", config);
    auto d = partitioner.partition(grid);

    EXPECT_EQ(d.nb_partitions(), 4);
    check_valid_partitions(d, 4);
    check_balanced(d, 4);
}

CASE("regional grid (non-global bounding box) partitions correctly") {
    // A regional lon-lat grid: the bounding box is a sub-region of the sphere,
    // exercising the general (non-global) bounding-box path.
    Grid grid(util::Config("type", "regional")("nx", 40)("ny", 30)("north", 60.)("south", 20.)("east", 40.)(
        "west", -20.));

    const int nb_parts = 6;
    grid::Partitioner partitioner("hilbert", nb_parts);
    auto d = partitioner.partition(grid);

    EXPECT_EQ(d.nb_partitions(), nb_parts);
    EXPECT_EQ(d.size(), static_cast<gidx_t>(grid.size()));
    check_valid_partitions(d, nb_parts);
    check_balanced(d, nb_parts);
}

CASE("space-filling curve preserves locality") {
    // Rough locality check: the mean squared distance between consecutive
    // grid points along the Hilbert ordering should be much smaller than
    // between random pairs. Here we verify a weaker but robust property:
    // each partition forms a spatially compact cluster, i.e. its bounding-box
    // area is significantly smaller than the whole-grid bounding-box area.
    Grid grid("O48");
    const int nb_parts = 8;
    grid::Partitioner partitioner("hilbert", nb_parts);
    auto d = partitioner.partition(grid);

    std::vector<double> xmin(nb_parts, 1e30), xmax(nb_parts, -1e30);
    std::vector<double> ymin(nb_parts, 1e30), ymax(nb_parts, -1e30);
    double gxmin = 1e30, gxmax = -1e30, gymin = 1e30, gymax = -1e30;

    gidx_t n = 0;
    for (const auto& ll : grid.lonlat()) {
        const int p = d.partition(n++);
        xmin[p]     = std::min(xmin[p], ll[0]);
        xmax[p]     = std::max(xmax[p], ll[0]);
        ymin[p]     = std::min(ymin[p], ll[1]);
        ymax[p]     = std::max(ymax[p], ll[1]);
        gxmin       = std::min(gxmin, ll[0]);
        gxmax       = std::max(gxmax, ll[0]);
        gymin       = std::min(gymin, ll[1]);
        gymax       = std::max(gymax, ll[1]);
    }

    const double global_area = (gxmax - gxmin) * (gymax - gymin);
    double sum_area          = 0.;
    for (int p = 0; p < nb_parts; ++p) {
        sum_area += (xmax[p] - xmin[p]) * (ymax[p] - ymin[p]);
    }
    // If partitions were spatially compact, the sum of their bounding-box areas
    // is at most a small multiple of the global area. A naive banded/round-robin
    // assignment would make each partition span the whole domain, giving
    // sum_area ~ nb_parts * global_area.
    EXPECT(sum_area < 2.0 * global_area);
}

//----------------------------------------------------------------------------------------------------------------------
// Weighted variant
//----------------------------------------------------------------------------------------------------------------------

CASE("unit weight function reproduces the unweighted partition (neutral default)") {
    Grid grid("O32");
    const int nb_parts = 4;

    SpaceFillingCurvePartitioner unweighted(nb_parts, util::NoConfig());
    SpaceFillingCurvePartitioner uniform(nb_parts, util::NoConfig());
    uniform.set_weight_function([](gidx_t, const PointXY&) { return 1.0; });

    EXPECT(not unweighted.has_weight_function());
    EXPECT(uniform.has_weight_function());

    std::vector<int> part_a(grid.size());
    std::vector<int> part_b(grid.size());
    unweighted.partition(grid, part_a.data());
    uniform.partition(grid, part_b.data());

    for (gidx_t n = 0; n < static_cast<gidx_t>(grid.size()); ++n) {
        EXPECT_EQ(part_a[n], part_b[n]);
    }
}

CASE("weighted partition balances total weight, not point count") {
    Grid grid("O48");
    const int nb_parts = 8;

    // A simple spatially-varying weight: northern-hemisphere points are 4x as
    // costly as southern-hemisphere points. This mimics an observation-density
    // imbalance concentrated in one region.
    const double heavy = 4.0;
    const double light = 1.0;
    auto weight        = [heavy, light](gidx_t, const PointXY& p) { return p.y() > 0. ? heavy : light; };

    SpaceFillingCurvePartitioner partitioner(nb_parts, util::NoConfig());
    partitioner.set_weight_function(weight);

    std::vector<int> part(grid.size());
    partitioner.partition(grid, part.data());

    // Accumulate per-partition weight and point count.
    std::vector<double> part_weight(nb_parts, 0.);
    std::vector<int> part_count(nb_parts, 0);
    double max_point_weight = 0.;
    gidx_t n                = 0;
    for (const auto& ll : grid.lonlat()) {
        const PointXY p{ll[0], ll[1]};
        const int pp = part[n++];
        EXPECT(pp >= 0);
        EXPECT(pp < nb_parts);
        const double w = (p.y() > 0. ? heavy : light);
        part_weight[pp] += w;
        part_count[pp] += 1;
        max_point_weight = std::max(max_point_weight, w);
    }

    const double w_min = *std::min_element(part_weight.begin(), part_weight.end());
    const double w_max = *std::max_element(part_weight.begin(), part_weight.end());
    const int c_min    = *std::min_element(part_count.begin(), part_count.end());
    const int c_max    = *std::max_element(part_count.begin(), part_count.end());

    // (a) Total weight is balanced: a greedy contiguous cut leaves each partition
    //     within one point's weight of the target on either side, so the spread is
    //     bounded by twice the heaviest point weight.
    EXPECT(w_max - w_min <= 2.0 * max_point_weight + 1.e-9);

    // (b) The weighting genuinely changed the balance away from equal point count:
    //     heavier (northern) regions get proportionally fewer points, so the point
    //     counts are now clearly unequal (unlike the unweighted case, where they
    //     differ by at most one).
    EXPECT(c_max > c_min + 1);
    EXPECT(c_max >= 2 * c_min);
}

CASE("masked (all-zero) weights fall back to a valid equal-count partition") {
    Grid grid("O32");
    const int nb_parts = 5;

    SpaceFillingCurvePartitioner partitioner(nb_parts, util::NoConfig());
    partitioner.set_weight_function([](gidx_t, const PointXY&) { return 0.0; });

    std::vector<int> part(grid.size());
    partitioner.partition(grid, part.data());

    std::vector<int> count(nb_parts, 0);
    for (gidx_t n = 0; n < static_cast<gidx_t>(grid.size()); ++n) {
        EXPECT(part[n] >= 0);
        EXPECT(part[n] < nb_parts);
        count[part[n]] += 1;
    }
    const int size       = static_cast<int>(grid.size());
    const int floor_size = size / nb_parts;
    const int ceil_size  = floor_size + (size % nb_parts != 0 ? 1 : 0);
    for (int c : count) {
        EXPECT(c >= floor_size);
        EXPECT(c <= ceil_size);
    }
}

CASE("weights can be supplied as a precomputed array (grid-order indexed)") {
    Grid grid("O48");
    const int nb_parts     = 8;
    const gidx_t nb_points = grid.size();

    // Build a precomputed weight array in grid-iteration order, equivalent to the
    // analytical northern-heavy weighting but expressed as raw per-point values -
    // e.g. an observation count per grid point.
    const double heavy = 4.0;
    const double light = 1.0;
    std::vector<double> weights(nb_points);
    {
        gidx_t n = 0;
        for (const auto& ll : grid.lonlat()) {
            weights[n++] = (ll[1] > 0. ? heavy : light);
        }
    }

    // Partition once via the analytical function and once via the array; the two
    // must produce identical partitions, demonstrating the array path is a faithful
    // alternative to the functional path.
    SpaceFillingCurvePartitioner by_function(nb_parts, util::NoConfig());
    by_function.set_weight_function([heavy, light](gidx_t, const PointXY& p) { return p.y() > 0. ? heavy : light; });

    SpaceFillingCurvePartitioner by_array(nb_parts, util::NoConfig());
    by_array.set_weights(weights);

    EXPECT(by_array.has_weight_function());

    std::vector<int> part_fn(nb_points);
    std::vector<int> part_arr(nb_points);
    by_function.partition(grid, part_fn.data());
    by_array.partition(grid, part_arr.data());

    for (gidx_t n = 0; n < nb_points; ++n) {
        EXPECT_EQ(part_fn[n], part_arr[n]);
    }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
