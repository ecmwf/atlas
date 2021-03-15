/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <mpi.h>

#include "eckit/exception/Exceptions.h"

namespace atlas {
  namespace repartition {
    namespace mpi {

      /// \brief Check if MPI routine succeeded.
      inline void callMpi(const int returnVal, const char* functionName,
        eckit::CodeLocation location) {

        // Check return value.
        if (returnVal != MPI_SUCCESS) {
          throw eckit::FailedLibraryCall("MPI " + std::string(MPICH_VERSION),
            functionName, "Error code: " + std::to_string(returnVal), location);
        }

        return;
      }
    #define CALL_MPI(a) callMpi(a, #a, Here())
    }
  }
}
