/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mpi/mpi.h"

namespace atlas {
namespace mpi {

extern "C"
{
  void atlas_mpi_comm_set_with_fortran_handle (int fhandle)
  {
    eckit::mpi::comm().set_with_fortran_handle(fhandle);
  }
  int atlas_mpi_comm_fortran_handle ()
  {
    return eckit::mpi::comm().fortran_handle();
  }
}

} // namespace mpi
} // namepsace atlas
