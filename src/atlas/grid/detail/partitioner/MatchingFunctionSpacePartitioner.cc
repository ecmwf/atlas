/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/MatchingFunctionSpacePartitioner.h"

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

MatchingFunctionSpacePartitioner::MatchingFunctionSpacePartitioner(): Partitioner() {
    ATLAS_NOTIMPLEMENTED;
}

// MatchingFunctionSpacePartitioner::MatchingFunctionSpacePartitioner(const idx_t nb_partitions):
//     Partitioner(nb_partitions) {
//     ATLAS_NOTIMPLEMENTED;
// }

MatchingFunctionSpacePartitioner::MatchingFunctionSpacePartitioner(const FunctionSpace& functionspace, const eckit::Parametrisation&):
    Partitioner(functionspace.nb_parts(),util::Config("mpi_comm",functionspace.mpi_comm())), partitioned_(functionspace) {}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
