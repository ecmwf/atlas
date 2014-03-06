#ifndef ATLAS_MPI_H
#define ATLAS_MPI_H

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif
 
#define MPI_COMM_WORLD 0

#define MPI_FAILURE 1
#define MPI_SUCCESS 0

#define MPI_STATUS_IGNORE 0
#define MPI_IN_PLACE ((void *) 1)
#define MPI_Comm int
#define MPI_Request int
#define MPI_Status int
#define MPI_Datatype int
#define MPI_Op int

#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_DOUBLE_PRECISION 3
#define MPI_BYTE 4

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
#define MPI_PRODUCT 4

void MPI_Abort ( MPI_Comm comm, int ierror );
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
int MPI_Bcast ( void *data, int n, MPI_Datatype datatype, int node, 
  MPI_Comm comm );
int MPI_Comm_rank ( MPI_Comm comm, int *me );
int MPI_Comm_size ( MPI_Comm comm, int *nprocs );
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
int MPI_Recv ( void *buf, int count, MPI_Datatype datatype,
  int source, int tag, MPI_Comm comm, MPI_Status *status );
int MPI_Reduce ( void *data1, void *data2, int n, MPI_Datatype datatype, 
  MPI_Op operation, int receiver, MPI_Comm comm );
int MPI_Reduce_scatter ( void *sendbuf, void *recvbuf, int recvcounts,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
int MPI_Rsend ( void *data, int n, MPI_Datatype datatype, int iproc, 
  int itag, MPI_Comm comm );
int MPI_Send ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm );
int MPI_Wait ( MPI_Request *request, MPI_Status *status );
int MPI_Waitall ( int icount, int irequest, MPI_Status status );
int MPI_Waitany ( int count, MPI_Request *request, int *index, 
  MPI_Status *status );

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif

#endif  // ATLAS_MPI_H