/*
 * (C) Copyright 1996-2013 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Tiago Quintino
/// @author Baudouin Raoult
/// @date   Nov 2013

#ifndef eckit_mpi_Comm_h
#define eckit_mpi_Comm_h

#include <mpi.h>

#include "eckit/exception/Exceptions.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace mpi {

//-----------------------------------------------------------------------------

/// MPI Exception
class Error: public eckit::Exception {
public:
    Error(const eckit::CodeLocation& loc, const std::string& s);
    Error(const eckit::CodeLocation& loc, const std::string& s, int errorcode);
};

//-----------------------------------------------------------------------------

/// Checking return values of MPI calls and throws exception on error
#define MPI_CHECK_RESULT( MPI_FUNC, args )                                                   \
{                                                                                            \
  int r = MPI_FUNC args;                                                                     \
  if (r != MPI_SUCCESS)                                                                      \
     throw eckit::mpi::Error(Here(),std::string(#MPI_FUNC#args),r);                           \
}

//-----------------------------------------------------------------------------

/// Singleton class to manage the MPI communication
class Comm {
public:

    /// Return a reference to the current PE
    static Comm& instance();

    /// @returns the generic communication channel
    MPI_Comm communicator();

    /// Initialise the PE
    /// @post will have a valid state
    void init(int argc=0, char** args=0);

    /// Free the PE, careful because some mpi implmentations fail upon re-init after a proper finalize
    /// @post will have not a valid state
    void finalize();

    /// @returns the MPI version as a pair of numbers
    std::pair<int,int> version() const;

    /// @returns the MPI version as a std::string
    std::string versionStr() const;

    /// Checks if the PE is initialized
    /// This is not the opposite of is_finalized
    bool isInitialised() const;

    /// Checks if the PE is finalized
    /// This is not the opposite of is_init
    bool isFinalised() const;

    /// Checks if the PE is in valid state
    /// Should be initialized and comm_ pointer is set
    bool isActive() const;

    /// Sets a barrier on a custom communicator.
    /// @param comm The communicator to set the barrier on.
    void barrier( MPI_Comm comm );

    /// Return rank, additionally, if is_init==0.
    size_t rank() const;

    /// Return the number of processes, or 1 if is_init==0.
    size_t size() const;

private: // methods

    Comm();

    ~Comm();

private:

    MPI_Comm comm_; ///< mpi communicator

};

//-----------------------------------------------------------------------------

} // namespace mpi
} // namespace eckit

#endif
