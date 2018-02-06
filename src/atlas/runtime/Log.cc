/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <chrono>
#include <thread>
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"
#include "eckit/os/BackTrace.h"

namespace atlas {

std::string backtrace() {
  return eckit::BackTrace::dump();
}

std::ostream& Log::debug_parallel() {
  const auto& comm = parallel::mpi::comm();
  auto rank = comm.rank();
  comm.barrier();
  for( int i=0; i<comm.size(); ++ i ) {
    std::this_thread::sleep_for(std::chrono::milliseconds( 10 ));
    if( i == rank ) {
      std::cout << "ATLAS_DEBUG_PARALLEL["<<rank<<"] ";
      return std::cout;
    }
  }
}

namespace detail {

void print_parallel_here(std::ostream& out, const eckit::CodeLocation& here) {
  const auto& comm = parallel::mpi::comm();
  comm.barrier();
  for( int i=0; i<comm.size(); ++ i ) {
    std::this_thread::sleep_for(std::chrono::milliseconds( 100 ));
    out << "[" << parallel::mpi::comm().rank() << "] DEBUG() @ " << here << std::endl;
  }
}
} // namespace detail

} // namespace atlas
