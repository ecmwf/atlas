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

#include "atlas/grid/detail/partitioner/MatchingFunctionSpacePartitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class MatchingFunctionSpacePartitionerLonLatPolygon : public MatchingFunctionSpacePartitioner {
public:
    static std::string static_type() { return "lonlat-polygon"; }

public:
    MatchingFunctionSpacePartitionerLonLatPolygon(): MatchingFunctionSpacePartitioner() {}
    // MatchingFunctionSpacePartitionerLonLatPolygon(const size_t nb_partitions):
    //     MatchingFunctionSpacePartitioner(nb_partitions) {}
    MatchingFunctionSpacePartitionerLonLatPolygon(const FunctionSpace& FunctionSpace, const eckit::Parametrisation& config):
        MatchingFunctionSpacePartitioner(FunctionSpace, config) {}

    using MatchingFunctionSpacePartitioner::partition;
    /**
   * @brief Partition a grid, using the same partitions from a pre-partitioned
   * FunctionSpace.
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
