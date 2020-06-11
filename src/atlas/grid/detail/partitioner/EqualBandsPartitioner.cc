/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "EqualBandsPartitioner.h"

#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

EqualBandsPartitioner::EqualBandsPartitioner() : EqualBandsPartitioner( mpi::size() ) {}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::EqualBandsPartitioner>
    __EqualBands( atlas::grid::detail::partitioner::EqualBandsPartitioner::static_type() );
}
