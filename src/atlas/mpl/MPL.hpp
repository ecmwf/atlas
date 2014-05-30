/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef MPL_hpp
#define MPL_hpp

#include <stdexcept>
#include <iostream>

#include "atlas/atlas_mpi.h"
#include "eckit/exception/Exceptions.h"

#define MPL_CHECK_RESULT( MPI_CALL )\
{ \
  int ierr = MPI_CALL; \
  if (ierr != MPI_SUCCESS) { \
    throw atlas::MPL::MPLError( std::string("MPI call: ") + \
      std::string(#MPI_CALL) + \
      std::string(" did not return MPI_SUCCESS."), Here() ); \
  } \
}

namespace atlas {

namespace MPL {


class MPLError : public eckit::Exception {
public:
  MPLError(const std::string& msg, const eckit::CodeLocation& loc)
  {
    eckit::StrStream s;
    s << "MPLError: " << msg << " " << " in " << loc << " "  << eckit::StrStream::ends;
    reason(std::string(s));
  }
};


  template<typename DATA_TYPE>
  MPI_Datatype TYPE();
  template<> inline MPI_Datatype TYPE<int>()    { return MPI_INT; }
  template<> inline MPI_Datatype TYPE<float>()  { return MPI_FLOAT; }
  template<> inline MPI_Datatype TYPE<double>() { return MPI_DOUBLE; }


  inline bool initialized()
  {
    int initialized;
    MPL_CHECK_RESULT( MPI_Initialized( &initialized ) );
    return initialized;
  }

  inline bool finalized()
  {
    int finalized;
    MPL_CHECK_RESULT( MPI_Finalized( &finalized ) );
    return finalized;
  }

  inline void init(int argc=0, char *argv[]=0)
  {
    if( !initialized() )
      MPL_CHECK_RESULT( MPI_Init(&argc,&argv) );
  }

  inline void finalize()
  {
    if( !finalized() )
      MPL_CHECK_RESULT( MPI_Finalize() );
  }

  inline int rank()
  {
    if( !initialized() ) throw MPLError( "MPI not initialized when calling MPL::rank()", Here() );
    int rank;
    MPL_CHECK_RESULT( MPI_Comm_rank( MPI_COMM_WORLD, &rank ) );
    return rank;
  }

  inline int size()
  {
    if( !initialized() ) throw MPLError( "MPI not initialized when calling MPL::size()", Here() );
    int nproc;
    MPL_CHECK_RESULT( MPI_Comm_size( MPI_COMM_WORLD, &nproc ) );
    return nproc;
  }

} // namespace MPL

} // namepsace atlas

#endif // MPL_hpp
