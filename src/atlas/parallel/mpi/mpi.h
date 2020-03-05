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

#include "eckit/mpi/Comm.h"

#include "atlas/parallel/mpi/Statistics.h"

namespace atlas {
namespace mpi {

using Comm = eckit::mpi::Comm;

inline const Comm& comm() {
    static const Comm& _comm = eckit::mpi::comm();
    return _comm;
}

inline idx_t rank() {
    return static_cast<idx_t>( comm().rank() );
}

inline int size() {
    return static_cast<idx_t>( comm().size() );
}

void finalize();
void finalise();

}  // namespace mpi
}  // namespace atlas
