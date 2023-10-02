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

#include "MatchingMeshPartitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class MatchingMeshPartitionerBruteForce : public MatchingMeshPartitioner {
public:
    static std::string static_type() { return "brute-force"; }

public:
    MatchingMeshPartitionerBruteForce(const Mesh& mesh, const eckit::Parametrisation& config): MatchingMeshPartitioner(mesh, config) {}


    MatchingMeshPartitionerBruteForce(): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerBruteForce(const eckit::Parametrisation&): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerBruteForce(const idx_t nb_partitions): MatchingMeshPartitioner() {}
    MatchingMeshPartitionerBruteForce(const idx_t nb_partitions, const eckit::Parametrisation&):
        MatchingMeshPartitioner() {}

    using MatchingMeshPartitioner::partition;
    virtual void partition(const Grid& grid, int partitioning[]) const;

    virtual std::string type() const { return static_type(); }
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
