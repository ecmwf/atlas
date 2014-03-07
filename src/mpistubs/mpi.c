#include "mpi.h"


#define MPI_INT_SIZE    sizeof(int)
#define MPI_FLOAT_SIZE  sizeof(float)
#define MPI_DOUBLE_SIZE sizeof(double)

int MPI_SIZE( MPI_Datatype datatype )
{
  switch( datatype ) {
    case(MPI_INT):    return MPI_INT_SIZE;
    case(MPI_FLOAT):  return MPI_FLOAT_SIZE;
    case(MPI_DOUBLE): return MPI_DOUBLE_SIZE;
  }
  return 0;
}

struct MPI_Request_t {
  MPI_Datatype datatype;
  int tag;
  int count;
  void *sendbuf;
  void *recvbuf;
} MPI_glb_request;



int MPI_Allgather ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype,
  MPI_Comm comm )
{
  memcpy( recvbuf, sendbuf, sendcount*MPI_SIZE(sendtype) );
  return MPI_SUCCESS;
}

int MPI_Allgatherv ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int *recvcounts, int *displs,
  MPI_Datatype recvtype, MPI_Comm comm )
{
  memcpy( recvbuf, sendbuf, sendcount*MPI_SIZE(sendtype) );
  return MPI_SUCCESS;
}

int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm ) 
{
  if( sendbuf == MPI_IN_PLACE ) return MPI_SUCCESS;
  memcpy( recvbuf, sendbuf, count*MPI_SIZE(datatype) );
  return MPI_SUCCESS;
}

int MPI_Alltoall (void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{ 
  memcpy( recvbuf, sendbuf, sendcount*MPI_SIZE(sendtype) );
  return MPI_SUCCESS; 
}

int MPI_Alltoallv (void *sendbuf, int *sendcounts, int *sdispls,
  MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
{ 
  memcpy( recvbuf, sendbuf, sendcounts[0]*MPI_SIZE(sendtype) );
  return MPI_SUCCESS; 
}

int MPI_Barrier ( MPI_Comm comm )
{ 
  return MPI_SUCCESS; 
}

int MPI_Comm_rank ( MPI_Comm comm, int *me ) 
{ 
  *me=0; 
  return MPI_SUCCESS; 
}

int MPI_Comm_size ( MPI_Comm comm, int *nprocs ) 
{ 
  *nprocs=1; 
  return MPI_SUCCESS; 
}

int MPI_Finalize ( void ) 
{ 
  return MPI_SUCCESS; 
}

int MPI_Finalized ( int *flag ) 
{ 
  *flag = 1; 
  return MPI_SUCCESS; 
}

int MPI_Gather (void *sendbuf, int sendcount, MPI_Datatype sendtype,
  void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) 
{ 
  memcpy( recvbuf, sendbuf, sendcount*MPI_SIZE(sendtype) );
  return MPI_SUCCESS; 
}

int MPI_Gatherv( void *sendbuf, int count, MPI_Datatype sendtype,
  void *recvbuf, int* recvcounts, int* displs, MPI_Datatype recvtype, int node,
  MPI_Comm comm )
{ 
  memcpy( recvbuf, sendbuf, count*MPI_SIZE(sendtype) );
  return MPI_SUCCESS; 
}

int MPI_Init ( int *argc, char **argv[] )
{ 
  return MPI_SUCCESS; 
}

int MPI_Initialized ( int *flag ) 
{ 
  *flag = 1; 
  return MPI_SUCCESS; 
}

int MPI_Irecv ( void *buf, int count, MPI_Datatype datatype,
  int source, int tag, MPI_Comm comm, MPI_Request *request )
{ 
  MPI_glb_request.recvbuf = buf;
  MPI_glb_request.count = count;
  MPI_glb_request.datatype = datatype;
  MPI_glb_request.tag = tag;
  return MPI_SUCCESS;
}

int MPI_Isend ( void *buf, int count, MPI_Datatype datatype,
  int dest, int tag, MPI_Comm comm, MPI_Request *request )
{ 
  MPI_glb_request.sendbuf = buf;
  MPI_glb_request.count = count;
  MPI_glb_request.datatype = datatype;
  MPI_glb_request.tag = tag;
  return MPI_SUCCESS;
}

int MPI_Wait ( MPI_Request *request, MPI_Status *status ) 
{ 
  memcpy( MPI_glb_request.recvbuf, MPI_glb_request.sendbuf, 
          MPI_glb_request.count*MPI_SIZE(MPI_glb_request.datatype) );
  if( status == MPI_STATUS_IGNORE ) return MPI_SUCCESS;
  *status=MPI_SUCCESS; 
  return MPI_SUCCESS; 
}
