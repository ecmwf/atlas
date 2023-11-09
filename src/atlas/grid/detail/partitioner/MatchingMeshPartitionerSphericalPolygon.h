/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class MatchingMeshPartitionerSphericalPolygon : public MatchingMeshPartitioner {
public:
    static std::string static_type() { return "spherical-polygon"; }

public:
    MatchingMeshPartitionerSphericalPolygon(const Mesh& mesh, const eckit::Parametrisation& config): MatchingMeshPartitioner(mesh, config) {}

    // Following throw ATLAS_NOT_IMPLEMENTED
    MatchingMeshPartitionerSphericalPolygon(): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerSphericalPolygon(const idx_t nb_partitions): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerSphericalPolygon(const idx_t nb_partitions, const eckit::Parametrisation& config):
        MatchingMeshPartitioner() {}
    MatchingMeshPartitionerSphericalPolygon(const eckit::Parametrisation&): MatchingMeshPartitioner() {}

    using MatchingMeshPartitioner::partition;
    /**
   * @brief Partition a grid, using the same partitions from a pre-partitioned
   * mesh.
   * The method constructs a partition edges polygon to test every target grid
   * node with.
   * @param[in] grid grid to be partitioned
   * @param[out] partitioning partitioning result
   */
    void partition(const Grid& grid, int partitioning[]) const;

    virtual std::string type() const { return static_type(); }
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
