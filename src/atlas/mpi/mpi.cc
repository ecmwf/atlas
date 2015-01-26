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

template<> MPI_Datatype datatype<char>()          { return MPI_CHAR; }
template<> MPI_Datatype datatype<int>()           { return MPI_INT; }
template<> MPI_Datatype datatype<unsigned int>()  { return MPI_UNSIGNED; }
template<> MPI_Datatype datatype<long>()          { return MPI_LONG; }
template<> MPI_Datatype datatype<unsigned long>() { return MPI_UNSIGNED_LONG; }
template<> MPI_Datatype datatype<float>()         { return MPI_FLOAT; }
template<> MPI_Datatype datatype<double>()        { return MPI_DOUBLE; }

template<> MPI_Datatype datatype<char>(char&)                   { return MPI_CHAR; }
template<> MPI_Datatype datatype<int>(int&)                     { return MPI_INT; }
template<> MPI_Datatype datatype<unsigned int>(unsigned int&)   { return MPI_UNSIGNED; }
template<> MPI_Datatype datatype<long>(long&)                   { return MPI_LONG; }
template<> MPI_Datatype datatype<unsigned long>(unsigned long&) { return MPI_UNSIGNED_LONG; }
template<> MPI_Datatype datatype<float>(float&)                 { return MPI_FLOAT; }
template<> MPI_Datatype datatype<double>(double&)               { return MPI_DOUBLE; }

bool initialized()
{
  int initialized;
  ATLAS_MPI_CHECK_RESULT( MPI_Initialized( &initialized ) );
  return initialized;
}

bool finalized()
{
  int finalized;
  ATLAS_MPI_CHECK_RESULT( MPI_Finalized( &finalized ) );
  return finalized;
}

void init(int argc, char *argv[])
{
  if( !initialized() )
    ATLAS_MPI_CHECK_RESULT( MPI_Init(&argc,&argv) );
}

void finalize()
{
  if( !finalized() )
    ATLAS_MPI_CHECK_RESULT( MPI_Finalize() );
}

Comm::Comm()
{
  if( ! initialized() )
    throw mpi::Error( "Trying to construct MPI communicator without MPI being initialized", Here() );
  set_with_C_handle(MPI_COMM_WORLD);
}

Comm::Comm( MPI_Comm comm )
{
  if( ! initialized() )
    throw mpi::Error( "Trying to construct MPI communicator without MPI being initialized", Here() );
  set_with_C_handle(comm);
}

Comm& Comm::instance()
{
  static Comm comm_instance;
  return comm_instance;
}

MPI_Fint Comm::fortran_handle()
{
  MPI_Fint fortran_comm = MPI_Comm_c2f(*this);
  return fortran_comm;
}

void Comm::set_with_fortran_handle( MPI_Fint handle )
{
  comm_ = MPI_Comm_f2c( handle );
}

void Comm::set_with_C_handle( MPI_Comm handle )
{
  comm_ = handle;
}

size_t Comm::rank() const
{
  if( !initialized() ) throw mpi::Error( "MPI not initialized when calling mpi::rank()", Here() );
  int rank;
  ATLAS_MPI_CHECK_RESULT( MPI_Comm_rank( *this, &rank ) );
  return rank;
}

size_t Comm::size() const
{
  if( !initialized() ) throw mpi::Error( "MPI not initialized when calling mpi::size()", Here() );
  int nproc;
  ATLAS_MPI_CHECK_RESULT( MPI_Comm_size( *this, &nproc ) );
  return nproc;
}

void Comm::barrier() const
{
  ATLAS_MPI_CHECK_RESULT( MPI_Barrier(*this) );
}

CommWorld& CommWorld::instance()
{
  static CommWorld world_instance;
  return world_instance;
}

size_t rank()
{
  return Comm::instance().rank();
}

size_t size()
{
  return Comm::instance().size();
}

void barrier()
{
  return Comm::instance().barrier();
}


extern "C"
{
  void atlas_mpi_comm_set_with_fortran_handle (int comm)
  {
    Comm::instance().set_with_fortran_handle(comm);
  }
  int atlas_mpi_comm_fortran_handle ()
  {
    return Comm::instance().fortran_handle();
  }
}

} // namespace mpi
} // namepsace atlas
