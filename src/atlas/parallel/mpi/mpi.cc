/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace mpi {

void finalize() {
    finalise();
}
void finalise() {
    Log::debug() << "atlas::mpi::finalize() --> Finalizing MPI" << std::endl;
    eckit::mpi::finaliseAllComms();
}

}  // namespace mpi
}  // namespace atlas
