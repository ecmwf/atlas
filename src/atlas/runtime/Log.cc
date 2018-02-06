/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"
#include "eckit/os/BackTrace.h"

namespace atlas {

std::string backtrace() {
  return eckit::BackTrace::dump();
}

namespace detail {

void debug_parallel_here( const eckit::CodeLocation& here ) {
  const auto& comm = parallel::mpi::comm();
  comm.barrier();
  Log::info() << "DEBUG_PARALLEL() @ " << here << std::endl;
}

void debug_parallel_what( const eckit::CodeLocation& here, const std::string& what ) {
  const auto& comm = parallel::mpi::comm();
  comm.barrier();
  Log::info() << "DEBUG_PARALLEL(" << what << ") @ " << here << std::endl;
}

} // namespace detail

} // namespace atlas
