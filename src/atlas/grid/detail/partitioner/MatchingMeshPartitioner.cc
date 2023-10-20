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

#include <sstream>

#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

MatchingMeshPartitioner::MatchingMeshPartitioner(): Partitioner() {
    ATLAS_THROW_EXCEPTION("Error: A MatchingMeshPartitioner needs to be initialised with a Mesh");
}

// MatchingMeshPartitioner::MatchingMeshPartitioner(const idx_t nb_partitions): Partitioner(nb_partitions) {
    // ATLAS_NOTIMPLEMENTED;
// }

MatchingMeshPartitioner::MatchingMeshPartitioner(const Mesh& mesh, const eckit::Parametrisation&):
    Partitioner(mesh.nb_parts(),util::Config("mpi_comm",mesh.mpi_comm())), prePartitionedMesh_(mesh) {}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
