/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/parallel/mpi/mpi.h"

#include <sstream>

#include "eckit/mpi/Comm.h"


namespace atlas {
namespace parallel {
namespace mpi {

const eckit::mpi::Comm& comm()
{
  return eckit::mpi::comm();
}

extern "C"
{
    void atlas_mpi_comm_attach_fortran_communicator (int fcomm )
    {
        std::ostringstream oss;
        oss << "fortran." << fcomm;
        eckit::mpi::addComm(oss.str().c_str(), fcomm );
    }

    int atlas_mpi_comm_fortran_communicator ()
    {
        NOTIMP;
    }
}

} // namespace mpi
} // namespace parallel
} // namespace atlas
