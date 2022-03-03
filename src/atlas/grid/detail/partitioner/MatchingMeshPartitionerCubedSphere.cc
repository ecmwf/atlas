/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid/Iterator.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.h"
#include "atlas/interpolation/method/cubedsphere/CellFinder.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

void MatchingMeshPartitionerCubedSphere::partition(const Grid& grid, int partitioning[]) const {

    // Make cell finder from owned mesh cells.
    // This is a little overkill, as it works out interpolation weights in the process.
    const auto finder = interpolation::method::cubedsphere::CellFinder(prePartitionedMesh_);

    // Loop over grid and set partioning[].
    auto lonlatIt = grid.lonlat().begin();
    for (gidx_t i = 0; i < grid.size(); ++i) {

        const auto lonlat = *lonlatIt;
        partitioning[i] = finder.getCell(lonlat).isect ? mpi::rank() : -1;;
        ++lonlatIt;
    }

    // AllReduce to get full partitioning array.
    mpi::comm().allReduceInPlace(partitioning, grid.size(), eckit::mpi::Operation::MAX);
    const auto misses = std::count_if(partitioning, partitioning + grid.size(), [](int elem){return elem < 0;});
    if (misses > 0) {
        throw_Exception(
            "Could not find partition for target node (source "
            "mesh does not contain all target grid points)\n" +
            std::to_string(misses) + " misses.", Here());
    }

}

} // namespace partitioner
} // namespace detail
} // namespace grid
} // namespace atlas
