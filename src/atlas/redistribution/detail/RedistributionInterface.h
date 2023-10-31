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

//#include "atlas/functionspace/detail/FunctionSpaceImpl.h"
#include "atlas/redistribution/Redistribution.h"

namespace atlas {
namespace functionspace{
class FunctionSpaceImpl;
}
namespace field {
class FieldSetImpl;
}
namespace grid {
class DistributionImpl;
}
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
class Grid;
}  // namespace grid
}  // namespace detail
}  // namespace grid
using GridImpl = grid::detail::grid::Grid;
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {
class Partitioner;
}  // namespace partitioner
}  // namespace detail
}  // namespace grid
using PartitionerImpl = grid::detail::partitioner::Partitioner;
}  // namespace atlas


namespace atlas {

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {

const Redistribution* atlas__Redistribution__new(
    const functionspace::FunctionSpaceImpl* fspace1, const functionspace::FunctionSpaceImpl* fspace2);

}

}  // namespace atlas
