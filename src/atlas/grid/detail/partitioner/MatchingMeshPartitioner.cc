/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

MatchingMeshPartitioner::MatchingMeshPartitioner(): Partitioner() {
    ATLAS_NOTIMPLEMENTED;
}

MatchingMeshPartitioner::MatchingMeshPartitioner(const idx_t nb_partitions): Partitioner(nb_partitions) {
    ATLAS_NOTIMPLEMENTED;
}

MatchingMeshPartitioner::MatchingMeshPartitioner(const Mesh& mesh):
    Partitioner(mesh.nb_partitions()), prePartitionedMesh_(mesh) {}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
