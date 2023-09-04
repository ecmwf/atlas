/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {


class MatchingMeshPartitionerCubedSphere : public MatchingMeshPartitioner {
public:
    static std::string static_type() { return "cubedsphere"; }

public:
    MatchingMeshPartitionerCubedSphere(): MatchingMeshPartitioner() {}
    // MatchingMeshPartitionerCubedSphere(const idx_t nb_partitions): MatchingMeshPartitioner(nb_partitions) {}
    // MatchingMeshPartitionerCubedSphere(const idx_t nb_partitions, const eckit::Parametrisation&):
    //     MatchingMeshPartitioner(nb_partitions) {}
    MatchingMeshPartitionerCubedSphere(const Mesh& mesh, const eckit::Parametrisation& config): MatchingMeshPartitioner(mesh, config) {}

    using MatchingMeshPartitioner::partition;
    virtual void partition(const Grid& grid, int partitioning[]) const;

    virtual std::string type() const { return static_type(); }
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
