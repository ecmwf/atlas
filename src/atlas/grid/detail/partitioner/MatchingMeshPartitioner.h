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

#include <vector>

#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class MatchingMeshPartitioner : public Partitioner {
public:
    MatchingMeshPartitioner();

    // MatchingMeshPartitioner(const idx_t nb_partitions);

    MatchingMeshPartitioner(const Mesh&, const eckit::Parametrisation&);

    virtual ~MatchingMeshPartitioner() override {}

protected:
    const Mesh prePartitionedMesh_;
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
