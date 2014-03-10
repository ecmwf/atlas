// (C) Copyright 1996-2014 ECMWF.

#ifndef MPL_hpp
#define MPL_hpp

#include <stdexcept>
#include <iostream>

#include "atlas/atlas_mpi.h"

#define MPL_CHECK_RESULT( MPI_CALL )\
{ \
  int ierr = MPI_CALL; \
  if (ierr != MPI_SUCCESS) { \
    throw std::runtime_error( std::string("MPI call: ") + \
      std::string(#MPI_CALL) + \
      std::string(" did not return MPI_SUCCESS.")); \
  } \
}

namespace atlas {

namespace MPL {

  template<typename DATA_TYPE>
  MPI_Datatype TYPE();
  template<> inline MPI_Datatype TYPE<int>()    { return MPI_INT; }
  template<> inline MPI_Datatype TYPE<float>()  { return MPI_FLOAT; }
  template<> inline MPI_Datatype TYPE<double>() { return MPI_DOUBLE; }

  

  inline void init(int argc=0, char *argv[]=0)
  {
    int initialized;
    MPL_CHECK_RESULT( MPI_Initialized( &initialized ) );
    if( !initialized )
      MPL_CHECK_RESULT( MPI_Init(&argc,&argv) );
  }

  inline void finalize()
  {
    int finalized;
    MPL_CHECK_RESULT( MPI_Finalized( &finalized ) );
    if( !finalized )
      MPL_CHECK_RESULT( MPI_Finalize() );
  }

  inline int rank()
  {
    int rank;
    MPL_CHECK_RESULT( MPI_Comm_rank( MPI_COMM_WORLD, &rank ) );
    return rank;
  }

  inline int size()
  {
    int nproc;
    MPL_CHECK_RESULT( MPI_Comm_size( MPI_COMM_WORLD, &nproc ) );
    return nproc;
  }

} // namespace MPL

} // namepsace atlas

#endif // MPL_hpp
