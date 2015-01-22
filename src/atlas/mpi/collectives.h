/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef ATLAS_MPI_COLLECTIVES_h
#define ATLAS_MPI_COLLECTIVES_h

#include "atlas/mpi/mpi.h"
#include "atlas/util/ArrayView.h"

namespace atlas {
namespace mpi {

/// @brief BufferBase
///
/// Class that keeps allocation of a MPI buffer including
/// counts and displacements
template <typename DATA_TYPE> struct BufferBase;
  
/// @brief Buffer<DATA_TYPE,SHAPE>
///
/// Class that keeps allocation of a MPI buffer including
/// counts and displacements, but with added index operator[]
/// that returns an ArrayView<DATA_TYPE,SHAPE> of the part
/// of the buffer for a processor index.
template <typename DATA_TYPE,int RANK> struct Buffer;

// ----------------------------------------------------------------------------------

template< typename DATA_TYPE >
inline int all_to_all( const Comm& comm,
                       std::vector< std::vector<DATA_TYPE> >& sendvec,
                       std::vector< std::vector<DATA_TYPE> >& recvvec );


template< typename DATA_TYPE >
inline int all_to_all( std::vector< std::vector<DATA_TYPE> >& sendvec,
                       std::vector< std::vector<DATA_TYPE> >& recvvec );


template<typename DATA_TYPE>
void all_gather( const Comm& comm, 
                 const DATA_TYPE send[], int sendcnt, 
                 BufferBase<DATA_TYPE>& recv );


template<typename DATA_TYPE>
void all_gather( const DATA_TYPE send[], int sendcnt, 
                 BufferBase<DATA_TYPE>& recv );


template<typename VECTOR>
void all_gather( const Comm& comm, 
                 const VECTOR& send, 
                 BufferBase<typename VECTOR::value_type>& recv );


template<typename VECTOR>
void all_gather( const VECTOR& send, 
                BufferBase<typename VECTOR::value_type>& recv );


// ----------------------------------------------------------------------------------

template <typename DATA_TYPE>
struct BufferBase
{
  int                    cnt;
  std::vector<int>       counts;
  std::vector<int>       displs;
  std::vector<DATA_TYPE> buf;

  BufferBase()
  {
    counts.resize( mpi::size() );
    displs.resize( mpi::size() );
  }
};

// ----------------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------------

template< typename DATA_TYPE >
inline int all_to_all( const Comm& comm,
                       std::vector< std::vector<DATA_TYPE> >& sendvec,
                       std::vector< std::vector<DATA_TYPE> >& recvvec )
{
  int cnt;
  int mpi_size = comm.size();
  // Get send-information
  std::vector<int> sendcounts(mpi_size);
  std::vector<int> senddispls(mpi_size);
  int sendcnt;
  senddispls[0] = 0;
  sendcounts[0] = sendvec[0].size();
  sendcnt = sendcounts[0];
  for( int jproc=1; jproc<mpi_size; ++jproc )
  {
    senddispls[jproc] = senddispls[jproc-1] + sendcounts[jproc-1];
    sendcounts[jproc] = sendvec[jproc].size();
    sendcnt += sendcounts[jproc];
  }


  // Get recv-information
  std::vector<int> recvcounts(mpi_size);
  std::vector<int> recvdispls(mpi_size);
  int recvcnt;
  ATLAS_MPI_CHECK_RESULT( MPI_Alltoall( sendcounts.data(), 1, MPI_INT,
                                  recvcounts.data(), 1, MPI_INT,
                                  comm ) );
  recvdispls[0] = 0;
  recvcnt = recvcounts[0];
  for( int jproc=1; jproc<mpi_size; ++jproc )
  {
    recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
    recvcnt += recvcounts[jproc];
  }

  // Communicate
  std::vector<DATA_TYPE> sendbuf(sendcnt);
  std::vector<DATA_TYPE> recvbuf(recvcnt);
  cnt = 0;
  for( int jproc=0; jproc<mpi_size; ++jproc )
  {
    for( int i=0; i<sendcounts[jproc]; ++i )
    {
      sendbuf[cnt++] = sendvec[jproc][i];
    }
  }
  ATLAS_MPI_CHECK_RESULT( MPI_Alltoallv(
                      sendbuf.data(), sendcounts.data(), senddispls.data(), mpi::datatype<DATA_TYPE>(),
                      recvbuf.data(), recvcounts.data(), recvdispls.data(), mpi::datatype<DATA_TYPE>(),
                      comm ) );
  cnt=0;
  for( int jproc=0; jproc<mpi_size; ++jproc )
  {
    recvvec[jproc].resize(recvcounts[jproc]);
    for( int i=0; i<recvcounts[jproc]; ++i )
    {
      recvvec[jproc][i] = recvbuf[cnt++];
    }
  }
  return MPI_SUCCESS;
}

template< typename DATA_TYPE >
inline int all_to_all( std::vector< std::vector<DATA_TYPE> >& sendvec,
                       std::vector< std::vector<DATA_TYPE> >& recvvec )
{
  return all_to_all( Comm::instance(), sendvec, recvvec );
}


template<typename DATA_TYPE>
inline void all_gather( const Comm& comm, const DATA_TYPE send[], int sendcnt, BufferBase<DATA_TYPE>& recv )
{
  int mpi_size = comm.size();
  ATLAS_MPI_CHECK_RESULT( MPI_Allgather( &sendcnt,           1, MPI_INT,
                                   recv.counts.data(), 1, MPI_INT, comm ) );
  recv.displs[0] = 0;
  recv.cnt = recv.counts[0];
  for( int jpart=1; jpart<mpi_size; ++jpart )
  {
    recv.displs[jpart] = recv.displs[jpart-1] + recv.counts[jpart-1];
    recv.cnt += recv.counts[jpart];
  }
  recv.buf.resize(recv.cnt);

  ATLAS_MPI_CHECK_RESULT( MPI_Allgatherv( const_cast<DATA_TYPE*>(send), sendcnt, mpi::datatype<DATA_TYPE>(),
                    recv.buf.data(), recv.counts.data(), recv.displs.data(),
                    mpi::datatype<DATA_TYPE>(), comm ) );
}
template<typename DATA_TYPE>
inline void all_gather( const DATA_TYPE send[], int sendcnt, BufferBase<DATA_TYPE>& recv )
{
  return all_gather( Comm::instance(), send, sendcnt, recv );
}

template<typename VECTOR>
inline void all_gather( const Comm& comm, const VECTOR& send, BufferBase<typename VECTOR::value_type>& recv )
{
  all_gather( comm, send.data(), send.size(), recv );
}

template<typename VECTOR>
inline void all_gather( const VECTOR& send, BufferBase<typename VECTOR::value_type>& recv )
{
  all_gather( Comm::instance(), send, recv );
}


} // namespace mpi
} // namepsace atlas

#endif // ATLAS_MPI_COLLECTIVES_h
