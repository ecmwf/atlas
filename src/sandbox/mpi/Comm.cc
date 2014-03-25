/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>

#include "eckit/log/Log.h"
#include "eckit/mpi/Comm.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace mpi {

static Mutex local_mutex_;

//-----------------------------------------------------------------------------

Error::Error(const eckit::CodeLocation& loc, const std::string& s)
    : Exception( s, loc )
{
}

Error::Error(const eckit::CodeLocation& loc, const std::string& s, int errorcode)
    : Exception( "", loc )
{
    std::ostringstream r;
    char errstr [MPI_MAX_ERROR_STRING];
    int errsize = 0;
    MPI_Error_string(errorcode,errstr,&errsize);
    std::string mpistr( errstr, errsize );

    r << "MPI Error: "
      << s
      << " ("
      << mpistr
      << " )" ;

    reason( r.str() ); // assign what
}

//-----------------------------------------------------------------------------

Comm::Comm()
{
    comm_ = MPI_COMM_NULL;
}

Comm::~Comm()
{
    finalize();
}

Comm& Comm::instance()
{
    AutoLock<Mutex> lock(local_mutex_);

    static Comm pe;
    return pe;
}

MPI_Comm Comm::communicator()
{
    ASSERT( isActive() );
    return comm_;
}

//-----------------------------------------------------------------------------

void Comm::init(int argc, char** args)
{
    AutoLock<Mutex> lock(local_mutex_);

    if ( isFinalised() )
        throw Error( Here(), "Should not call Comm::initialize() after Comm::finalize()" );

    if( !isInitialised() ) // then initialize
    {
        MPI_CHECK_RESULT(MPI_Init,(&argc,&args));
    }

    comm_ = MPI_COMM_WORLD;
}

//-----------------------------------------------------------------------------

void Comm::finalize()
{
    AutoLock<Mutex> lock(local_mutex_);

    if( isInitialised() && !isFinalised() )
    {
        MPI_CHECK_RESULT(MPI_Finalize,());
    }

    comm_ = MPI_COMM_NULL;
}

//-----------------------------------------------------------------------------

bool Comm::isInitialised() const
{
    AutoLock<Mutex> lock(local_mutex_);

    int is_init = 0;
    MPI_CHECK_RESULT(MPI_Initialized,(&is_init));
    return bool(is_init);
}

//-----------------------------------------------------------------------------

bool Comm::isFinalised() const
{
    AutoLock<Mutex> lock(local_mutex_);

    int is_finalized = 0;
    MPI_CHECK_RESULT(MPI_Finalized,(&is_finalized));
    return bool(is_finalized);
}

//-----------------------------------------------------------------------------

bool Comm::isActive() const
{
    AutoLock<Mutex> lock(local_mutex_);

    return isInitialised() && !isFinalised()  && ( comm_ != MPI_COMM_NULL );
}

//-----------------------------------------------------------------------------

std::pair<int,int> Comm::version() const
{
    int version = 0;
    int subversion = 0;
    MPI_CHECK_RESULT(MPI_Get_version,(&version,&subversion));
    return std::make_pair( version, subversion );
}

//-----------------------------------------------------------------------------

std::string Comm::versionStr() const
{
    std::pair<int,int> v = version();
    std::ostringstream s;
    s << v.first
      << "."
      << v.second;
    return s.str();
}

//-----------------------------------------------------------------------------

void Comm::barrier( MPI_Comm comm )
{
    AutoLock<Mutex> lock(local_mutex_);

    ASSERT( comm != MPI_COMM_NULL );
    ASSERT( isActive() );

    MPI_CHECK_RESULT(MPI_Barrier,(comm));
}

//-----------------------------------------------------------------------------

size_t Comm::rank() const
{
    AutoLock<Mutex> lock(local_mutex_);

    if ( !isActive() ) return 0;
    int irank;
    MPI_CHECK_RESULT(MPI_Comm_rank,(comm_,&irank));
    return static_cast<size_t>(irank);
}

//-----------------------------------------------------------------------------

size_t Comm::size() const
{
    AutoLock<Mutex> lock(local_mutex_);

    if ( !isActive() ) return 1;
    int nproc;
    MPI_CHECK_RESULT(MPI_Comm_size,(comm_,&nproc));
    return static_cast<size_t>(nproc);
}

//-----------------------------------------------------------------------------

} // namespace mpi
} // namespace eckit
