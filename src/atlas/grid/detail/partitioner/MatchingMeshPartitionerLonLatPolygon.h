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

class MatchingMeshPartitionerLonLatPolygon : public MatchingMeshPartitioner {
public:
    static std::string static_type() { return "lonlat-polygon"; }

public:
    MatchingMeshPartitionerLonLatPolygon(const Mesh& mesh, const eckit::Parametrisation& config);


    MatchingMeshPartitionerLonLatPolygon(): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerLonLatPolygon(const eckit::Parametrisation&): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerLonLatPolygon(const size_t nb_partitions): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerLonLatPolygon(const size_t nb_partitions, const eckit::Parametrisation& config):
        MatchingMeshPartitioner() {}

    using Partitioner::partition;

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

private:
    bool fallback_nearest_{false};
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
