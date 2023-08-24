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

#include <string_view>

#include "eckit/mpi/Comm.h"
#include "atlas/parallel/mpi/Statistics.h"

namespace atlas {
namespace mpi {

using Comm = eckit::mpi::Comm;

inline const Comm& comm() {
    return eckit::mpi::comm();
}

inline const Comm& comm(std::string_view name) {
    return eckit::mpi::comm(name.data());
}

inline idx_t rank() {
    return static_cast<idx_t>(comm().rank());
}

inline int size() {
    return static_cast<idx_t>(comm().size());
}

void finalize();
void finalise();

}  // namespace mpi
}  // namespace atlas
