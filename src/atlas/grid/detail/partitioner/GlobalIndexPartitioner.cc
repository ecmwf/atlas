/*
 * (C) Crown copyright 2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "GlobalIndexPartitioner.h"

#include <vector>

#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

GlobalIndexPartitioner::GlobalIndexPartitioner():
    Partitioner() {
}

GlobalIndexPartitioner::GlobalIndexPartitioner(int N):
    GlobalIndexPartitioner(N, util::NoConfig()) {
}

GlobalIndexPartitioner::GlobalIndexPartitioner(int N, const eckit::Parametrisation& config):
    Partitioner(N, config) {
}

GlobalIndexPartitioner::GlobalIndexPartitioner(const eckit::Parametrisation& config):
    Partitioner(config) {
}

void GlobalIndexPartitioner::partition(const Grid& grid, int part[]) const {
    const size_t quotient = grid.size() / atlas::mpi::comm().size();
    const size_t remainder = grid.size() % atlas::mpi::comm().size();

    for (auto n = 0; n < grid.size(); n++) {
        // PEs 0 to (remainder - 1) have (quotient + 1) points
        // Remaining PEs have (quotient) points
        auto r = n / (quotient + 1);
        part[n] = (r >= remainder) ? (n - remainder) / quotient : r;
    }
}

Distribution GlobalIndexPartitioner::partition(const Grid& grid) const {
    std::vector<int> partitions;
    partitions.reserve(grid.size());

    const size_t quotient = grid.size() / atlas::mpi::comm().size();
    const size_t remainder = grid.size() % atlas::mpi::comm().size();

    for (auto n = 0; n < grid.size(); n++) {
        // PEs 0 to (remainder - 1) have (quotient + 1) points
        // Remaining PEs have (quotient) points
        auto r = n / (quotient + 1);
        partitions[n] = (r >= remainder) ? (n - remainder) / quotient : r;
    }

    return Distribution(atlas::mpi::comm().size(), grid.size(), partitions.data());
}

}
}
}
}

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::GlobalIndexPartitioner>
    __GlobalIndex("global_index");
}
