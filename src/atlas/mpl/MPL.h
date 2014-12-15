/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef MPL_h
#define MPL_h

#include <stdexcept>
#include <iostream>

#include "atlas/atlas_mpi.h"
#include "atlas/util/ArrayView.h"
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
  template<> inline MPI_Datatype TYPE<char>()          { return MPI_CHAR; }
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


  template <typename DATA_TYPE>
  struct BufferBase
  {
    int                    cnt;
    std::vector<int>       counts;
    std::vector<int>       displs;
    std::vector<DATA_TYPE> buf;

    BufferBase()
    {
      counts.resize( MPL::size() );
      displs.resize( MPL::size() );
    }
  };

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


  template<typename DATA_TYPE>
  void all_gather( const DATA_TYPE send[], int sendcnt, BufferBase<DATA_TYPE>& recv )
  {
    MPL_CHECK_RESULT( MPI_Allgather( &sendcnt,           1, MPI_INT,
                                     recv.counts.data(), 1, MPI_INT, MPI_COMM_WORLD ) );
    recv.displs[0] = 0;
    recv.cnt = recv.counts[0];
    for( int jpart=1; jpart<MPL::size(); ++jpart )
    {
      recv.displs[jpart] = recv.displs[jpart-1] + recv.counts[jpart-1];
      recv.cnt += recv.counts[jpart];
    }
    recv.buf.resize(recv.cnt);

    MPL_CHECK_RESULT( MPI_Allgatherv( const_cast<DATA_TYPE*>(send), sendcnt, MPL::TYPE<DATA_TYPE>(),
                      recv.buf.data(), recv.counts.data(), recv.displs.data(),
                      MPL::TYPE<DATA_TYPE>(), MPI_COMM_WORLD) );
  }

  template<typename VECTOR>
  void all_gather( const VECTOR& send, BufferBase<typename VECTOR::value_type>& recv )
  {
    all_gather(send.data(),send.size(),recv);
  }


  template <typename DATA_TYPE,int RANK>
  struct Buffer : BufferBase<DATA_TYPE>
  {
  };

  template <typename DATA_TYPE>
  struct Buffer<DATA_TYPE,1> : public BufferBase<DATA_TYPE>
  {
    ArrayView<DATA_TYPE,1> operator[](int p)
    {
      return ArrayView<DATA_TYPE,1> ( BufferBase<DATA_TYPE>::buf.data()+BufferBase<DATA_TYPE>::displs[p],
                                      make_shape( BufferBase<DATA_TYPE>::counts[p] ).data() );
    }
  };


} // namespace MPL

} // namepsace atlas

#endif // MPL_h
