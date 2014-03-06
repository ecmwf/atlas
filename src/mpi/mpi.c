#include "mpi.h"

int MPI_Barrier ( MPI_Comm comm ) { return MPI_SUCCESS; }
int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm ) { return MPI_FAILURE; }
int MPI_Alltoall (void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) { return MPI_FAILURE; }
int MPI_Alltoallv (void *sendbuf, int *sendcounts, int *sdispls,
  MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm) { return MPI_FAILURE; }
int MPI_Comm_rank ( MPI_Comm comm, int *me ) { *me=0; return MPI_SUCCESS; }
int MPI_Comm_size ( MPI_Comm comm, int *nprocs ) { *nprocs=1; return MPI_SUCCESS; }
int MPI_Gather (void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) { return MPI_FAILURE; }
int MPI_Init ( int *argc, char **argv[] ) { return MPI_SUCCESS; }
int MPI_Initialized ( int *flag ) { *flag = 1; return MPI_SUCCESS; }
int MPI_Finalize ( void ) { return MPI_SUCCESS; }
int MPI_Finalized ( int *flag ) { *flag = 1; return MPI_SUCCESS; }
int MPI_Irecv ( void *buf, int count, MPI_Datatype datatype,
  int source, int tag, MPI_Comm comm, MPI_Request *request ) { return MPI_FAILURE; }
int MPI_Isend ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm, MPI_Request *request ) { return MPI_FAILURE; }
int MPI_Gatherv( void *sendbuf, int count, MPI_Datatype senddatatype,
  void *recvbuf, int* recvcounts, int* displs, MPI_Datatype recvdatatype, int node,
  MPI_Comm comm ) { return MPI_FAILURE; }
int MPI_Wait ( MPI_Request *request, MPI_Status *status ) { *status=MPI_SUCCESS; return MPI_SUCCESS; }
