/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation/method/cubedsphere/CellFinder.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

void MatchingMeshPartitionerCubedSphere::partition(const Grid& grid, int partitioning[]) const {
    // Make cell finder from owned mesh cells.
    const auto finder = interpolation::method::cubedsphere::CellFinder(prePartitionedMesh_);

    // Numeric tolerance should scale with N.
    const auto N           = CubedSphereGrid(prePartitionedMesh_.grid()).N();
    const auto epsilon     = 2. * std::numeric_limits<double>::epsilon() * N;
    const auto edgeEpsilon = epsilon;
    const size_t listSize  = 8;

    // Loop over grid and set partioning[].
    auto lonlatIt = grid.lonlat().begin();
    for (gidx_t i = 0; i < grid.size(); ++i) {
        // This is probably more expensive than it needs to be, as it performs
        // a dry run of the cubedsphere interpolation method.
        const auto& lonlat = *lonlatIt;
        partitioning[i]    = finder.getCell(lonlat, listSize, edgeEpsilon, epsilon).isect ? mpi::rank() : -1;
        ++lonlatIt;
    }

    // AllReduce to get full partitioning array.
    mpi::comm().allReduceInPlace(partitioning, grid.size(), eckit::mpi::Operation::MAX);
    const auto misses = std::count_if(partitioning, partitioning + grid.size(), [](int elem) { return elem < 0; });
    if (misses > 0) {
        throw_Exception(
            "Could not find partition for target node (source "
            "mesh does not contain all target grid points)\n" +
                std::to_string(misses) + " misses.",
            Here());
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
