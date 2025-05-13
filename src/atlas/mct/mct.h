/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/mpi/Comm.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas::mct {

  void set_default_comm_to_local(int model_id) {
    mpi::comm().split(model_id, std::to_string(model_id));
    eckit::mpi::setCommDefault(std::to_string(model_id));
  }

} // namespace atlas::mct
