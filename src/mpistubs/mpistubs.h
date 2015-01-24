/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef MPISTUBS_H
#define MPISTUBS_H

/* If this is a C++ compiler, use C linkage */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif
 
#define MPI_COMM_WORLD 0

#define MPI_FAILURE 1
#define MPI_SUCCESS 0
#define MPI_MAX_ERROR_STRING 0

#define MPI_Status int
#define MPI_STATUS_IGNORE ((MPI_Status *) 0)
#define MPI_IN_PLACE ((void *) 1)
typedef struct mpi_communicator_t { int comm; } *MPI_Comm;

#define MPI_Datatype int
#define MPI_Op int

#define MPI_CHAR 0
#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_DOUBLE_PRECISION 3
#define MPI_BYTE 4
#define MPI_UNSIGNED 5
#define MPI_UNSIGNED_LONG 6
#define MPI_LONG 7
#define MPI_Fint int

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
#define MPI_PRODUCT 4

#define MPI_Request int

int MPI_Comm_test_inter ( MPI_Comm comm, int *flag );
int MPI_Allgather ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype,
  MPI_Comm comm );
int MPI_Allgatherv ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int *recvcounts, int *displs,
  MPI_Datatype recvtype, MPI_Comm comm );
int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
  MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Barrier ( MPI_Comm comm );
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
  int root, MPI_Comm comm);
MPI_Fint MPI_Comm_c2f( MPI_Comm comm );
MPI_Comm MPI_Comm_f2c( MPI_Fint f_handle );
int MPI_Comm_rank ( MPI_Comm comm, int *me );
int MPI_Comm_size ( MPI_Comm comm, int *nprocs );
int MPI_Error_string ( int errcode, char* errstr, int *errsize );
int MPI_Finalize ( void );
int MPI_Finalized ( int *flag );
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Gatherv( void *sendbuf, int count, MPI_Datatype senddatatype,
  void *recvbuf, int* recvcounts, int* displs, MPI_Datatype recvdatatype, 
  int node, MPI_Comm comm );
int MPI_Init ( int *argc, char **argv[] );
int MPI_Initialized ( int *flag );
int MPI_Irecv ( void *buf, int count, MPI_Datatype datatype,
  int source, int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Isend ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm, MPI_Request *request );
int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, MPI_Comm comm);
int MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Wait ( MPI_Request *request, MPI_Status *status );

/* If this is a C++ compiler, end C linkage */
#if defined(c_plusplus) || defined(__cplusplus) 
}
#endif

#endif  // MPISTUBS_H
