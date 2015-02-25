/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef ATLAS_MPI_MPI_h
#define ATLAS_MPI_MPI_h

#include "eckit/mpi/mpi.h"
#include "eckit/mpi/Collectives.h"
#include "eckit/mpi/Exceptions.h"
#include "atlas/mpi/Collectives.h"

namespace atlas {
namespace mpi {

} // namespace mpi
} // namespace atlas

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas_mpi_comm_attach_fortran_communicator (int comm);
  int atlas_mpi_comm_fortran_communicator ();
}
// ------------------------------------------------------------------

#endif // ATLAS_MPI_MPI_h
