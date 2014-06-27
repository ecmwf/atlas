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
  template<> inline MPI_Datatype TYPE<int>()           { return MPI_INT; }
  template<> inline MPI_Datatype TYPE<unsigned int>()  { return MPI_UNSIGNED; }
  template<> inline MPI_Datatype TYPE<long>()          { return MPI_LONG; }
  template<> inline MPI_Datatype TYPE<unsigned long>() { return MPI_UNSIGNED_LONG; }
  template<> inline MPI_Datatype TYPE<float>()         { return MPI_FLOAT; }
  template<> inline MPI_Datatype TYPE<double>()        { return MPI_DOUBLE; }


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

  template< typename DATA_TYPE >
  inline int Alltoall( std::vector< std::vector<DATA_TYPE> >& sendvec,
                       std::vector< std::vector<DATA_TYPE> >& recvvec )
  {
    int cnt;

    // Get send-information
    std::vector<int> sendcounts(MPL::size());
    std::vector<int> senddispls(MPL::size());
    int sendcnt;
    senddispls[0] = 0;
    sendcounts[0] = sendvec[0].size();
    sendcnt = sendcounts[0];
    for( int jproc=1; jproc<MPL::size(); ++jproc )
    {
      senddispls[jproc] = senddispls[jproc-1] + sendcounts[jproc-1];
      sendcounts[jproc] = sendvec[jproc].size();
      sendcnt += sendcounts[jproc];
    }


    // Get recv-information
    std::vector<int> recvcounts(MPL::size());
    std::vector<int> recvdispls(MPL::size());
    int recvcnt;
    MPL_CHECK_RESULT( MPI_Alltoall( sendcounts.data(), 1, MPI_INT,
                                    recvcounts.data(), 1, MPI_INT,
                                    MPI_COMM_WORLD ) );
    recvdispls[0] = 0;
    recvcnt = recvcounts[0];
    for( int jproc=1; jproc<MPL::size(); ++jproc )
    {
      recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
      recvcnt += recvcounts[jproc];
    }

    // Communicate
    std::vector<DATA_TYPE> sendbuf(sendcnt);
    std::vector<DATA_TYPE> recvbuf(recvcnt);
    cnt = 0;
    for( int jproc=0; jproc<MPL::size(); ++jproc )
    {
      for( int i=0; i<sendcounts[jproc]; ++i )
      {
        sendbuf[cnt++] = sendvec[jproc][i];
      }
    }
    MPL_CHECK_RESULT( MPI_Alltoallv(
                        sendbuf.data(), sendcounts.data(), senddispls.data(), MPL::TYPE<DATA_TYPE>(),
                        recvbuf.data(), recvcounts.data(), recvdispls.data(), MPL::TYPE<DATA_TYPE>(),
                        MPI_COMM_WORLD ) );
    cnt=0;
    for( int jproc=0; jproc<MPL::size(); ++jproc )
    {
      recvvec[jproc].resize(recvcounts[jproc]);
      for( int i=0; i<recvcounts[jproc]; ++i )
      {
        recvvec[jproc][i] = recvbuf[cnt++];
      }
    }
    return MPI_SUCCESS;
  }

} // namespace MPL

} // namepsace atlas

#endif // MPL_hpp
