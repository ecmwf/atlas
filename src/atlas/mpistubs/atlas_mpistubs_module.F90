
module atlas_mpistubs_module
  use, intrinsic :: iso_c_binding, only : c_float, c_double
  implicit none
  public
!
!  Dummy parameters for MPI F90 stubs
!
  integer mpi_comm_world
  parameter ( mpi_comm_world = 0 )
!
!  Return values
!
  integer mpi_failure
  parameter ( mpi_failure = 1 )
  integer mpi_success
  parameter ( mpi_success = 0 )
!
!  recv message status
!
  integer mpi_status_size
  parameter ( mpi_status_size = 4 )
  integer mpi_source
  parameter ( mpi_source = 1 )
  integer mpi_tag
  parameter ( mpi_tag = 2 )
  integer mpi_count
  parameter ( mpi_count = 3 )
!
!  recv flags
!
  integer mpi_any_source
  parameter ( mpi_any_source = -1 )
  integer mpi_any_tag
  parameter ( mpi_any_tag = -1 )
!
!  data types and sizes
!
  integer MPI_INTEGER
  parameter ( MPI_INTEGER = 28 )
  integer MPI_REAL
  parameter ( MPI_REAL = 11 )
  integer MPI_REAL4
  parameter ( MPI_REAL4 = MPI_REAL )
  integer mpi_double_precision
  parameter ( mpi_double_precision = 3 )
  integer MPI_REAL8
  parameter ( MPI_REAL8 = mpi_double_precision )
  integer mpi_logical
  parameter ( mpi_logical = 4 )
  integer mpi_character
  parameter ( mpi_character = 5 )
  integer MPI_IN_PLACE
  parameter ( MPI_IN_PLACE = 0 )

!
!  allreduce operations
!
  INTEGER MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD
  PARAMETER (MPI_MAX=100,MPI_MIN=101,MPI_SUM=102,MPI_PROD=103)

  interface mpi_gather
     module procedure mpi_gather_int_r0
     module procedure mpi_gather_int_r1
     module procedure mpi_gather_real32_r0
     module procedure mpi_gather_real32_r1
     module procedure mpi_gather_real64_r0
     module procedure mpi_gather_real64_r1
  end interface mpi_gather

  interface mpi_gatherv
     module procedure mpi_gatherv_int_r1
     module procedure mpi_gatherv_real32_r1
     module procedure mpi_gatherv_real64_r1
  end interface mpi_gatherv

  interface mpi_allgather
     module procedure mpi_allgather_int_r1
     module procedure mpi_allgather_real32_r1
     module procedure mpi_allgather_real64_r1
  end interface mpi_allgather
  
  interface mpi_allgatherv
     module procedure mpi_allgatherv_int_r1
     module procedure mpi_allgatherv_real32_r1
     module procedure mpi_allgatherv_real64_r1
  end interface mpi_allgatherv
  
  interface mpi_alltoall
    module procedure mpi_alltoall_int_r1
    module procedure mpi_alltoall_real32_r1
    module procedure mpi_alltoall_real64_r1
  end interface mpi_alltoall
  
  interface mpi_alltoallv
    module procedure mpi_alltoallv_int_r1
    module procedure mpi_alltoallv_real32_r1
    module procedure mpi_alltoallv_real64_r1
  end interface mpi_alltoallv
    
  interface mpi_allreduce
     module procedure mpi_allreduce_int_r0
     module procedure mpi_allreduce_int_r1
     module procedure mpi_allreduce_real32_r0
     module procedure mpi_allreduce_real32_r1
     module procedure mpi_allreduce_real64_r0
     module procedure mpi_allreduce_real64_r1
     module procedure mpi_allreduce_int_r1_inplace
     module procedure mpi_allreduce_real32_r0_inplace
     module procedure mpi_allreduce_real32_r1_inplace
  end interface mpi_allreduce
  
  interface mpi_reduce
     module procedure mpi_reduce_int_r0
     module procedure mpi_reduce_int_r1
     module procedure mpi_reduce_real32_r0
     module procedure mpi_reduce_real64_r0
     module procedure mpi_reduce_real32_r1
     module procedure mpi_reduce_real64_r1
  end interface mpi_reduce
  
  interface mpi_reduce_scatter
     module procedure mpi_reduce_scatter_int_r1
     module procedure mpi_reduce_scatter_real32_r1
  end interface mpi_reduce_scatter
  
  interface mpi_bcast
     module procedure mpi_bcast_int_r0
     module procedure mpi_bcast_int_r1
     module procedure mpi_bcast_int_r2
     module procedure mpi_bcast_real32_r0
     module procedure mpi_bcast_real64_r0
     module procedure mpi_bcast_real32_r1
     module procedure mpi_bcast_real64_r1
     module procedure mpi_bcast_real32_r2
     module procedure mpi_bcast_real64_r2
     module procedure mpi_bcast_char
  end interface mpi_bcast

  interface mpi_isend
     module procedure mpi_isend_int_r0
     module procedure mpi_isend_int_r1
     module procedure mpi_isend_real32_r0
     module procedure mpi_isend_real32_r1
     module procedure mpi_isend_real64_r0
     module procedure mpi_isend_real64_r1
  end interface mpi_isend

  interface mpi_irecv
     module procedure mpi_irecv_int_r0
     module procedure mpi_irecv_int_r1
     module procedure mpi_irecv_real32_r0
     module procedure mpi_irecv_real32_r1
     module procedure mpi_irecv_real64_r0
     module procedure mpi_irecv_real64_r1
  end interface mpi_irecv
  
  interface mpi_send
     module procedure mpi_send_int_r0
     module procedure mpi_send_int_r1
     module procedure mpi_send_int_r2
     module procedure mpi_send_int_r2_fixedstatus
     module procedure mpi_send_real32_r0
     module procedure mpi_send_real32_r1
     module procedure mpi_send_real64_r0
     module procedure mpi_send_real64_r1
  end interface mpi_send
  
  interface mpi_recv
     module procedure mpi_recv_int_r0
     module procedure mpi_recv_int_r0_fixedstatus
     module procedure mpi_recv_int_r1
     module procedure mpi_recv_int_r2
     module procedure mpi_recv_real32_r0
     module procedure mpi_recv_real32_r1
  end interface mpi_recv
  
  interface mpi_pack
     module procedure mpi_pack_int_r0
     module procedure mpi_pack_int_r1
     module procedure mpi_pack_real32_r0
     module procedure mpi_pack_real32_r1
  end interface mpi_pack

  interface mpi_unpack
     module procedure mpi_unpack_int_r0
     module procedure mpi_unpack_int_r1
     module procedure mpi_unpack_real32_r0
     module procedure mpi_unpack_real32_r1
  end interface mpi_unpack
  
  interface mpi_copy
     module procedure mpi_copy_int_r1
     module procedure mpi_copy_real32_r1
     module procedure mpi_copy_real64_r1
  end interface mpi_copy
  

CONTAINS
  
  !*****************************************************************************
  !
  ! MPI_ABORT shuts down the processes in a given communicator.
  !
  !*****************************************************************************
  
  subroutine mpi_abort ( comm, errorcode, ierror )
    
    implicit none
    
    integer comm, errorcode, ierror
    integer, parameter :: MPI_FAILURE = 1
    integer, parameter :: MPI_SUCCESS = 0
    
    !ierror = MPI_SUCCESS
    
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_ABORT:'
    write ( *, '(a,i12)' ) '  Shut down with error code = ', errorcode
    
    stop
  end subroutine mpi_abort
  
  !*****************************************************************************
  !
  ! MPI_GATHER gathers data from all the processes in a communicator.
  !
  !*****************************************************************************
  
  subroutine mpi_gather_int_r0( data1, nsend, sendtype, data2, nrecv, recvtype, root, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror, root
    integer :: data1, data2(nsend)
    
    ierror = MPI_SUCCESS
    
    data2(:) = data1
    
  end subroutine mpi_gather_int_r0
  
  subroutine mpi_gather_real32_r0( data1, nsend, sendtype, data2, nrecv, recvtype, root, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror, root
    real :: data1, data2(nsend)
    
    ierror = MPI_SUCCESS
    
    data2(:) = data1
    
  end subroutine mpi_gather_real32_r0
  
  subroutine mpi_gather_real64_r0( data1, nsend, sendtype, data2, nrecv, recvtype, root, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror, root
    real(c_double) :: data1, data2(nsend)
    
    ierror = MPI_SUCCESS
    data2(:) = data1
    
  end subroutine mpi_gather_real64_r0
  
  subroutine mpi_gather_int_r1( data1, nsend, sendtype, data2, nrecv, recvtype, root, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror, root
    integer :: data1(nsend), data2(nsend)
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_gather_int_r1
  
  subroutine mpi_gather_real32_r1( data1, nsend, sendtype, data2, nrecv, recvtype, root, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror, root
    real :: data1(nsend), data2(nsend)
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_gather_real32_r1
  
  subroutine mpi_gather_real64_r1( data1, nsend, sendtype,data2, nrecv, recvtype, root, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror, root
    real(c_double) :: data1(nsend), data2(nsend)
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_gather_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_ALLGATHER gathers data from all the processes in a communicator.
  !
  !*****************************************************************************
  
  subroutine mpi_allgather_int_r1( data1, nsend, sendtype, data2, nrecv, recvtype, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror
    integer :: data1(nsend), data2(nsend)
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_allgather_int_r1
  
  subroutine mpi_allgather_real32_r1( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror
    real(c_float) :: data1(nsend), data2(nsend)
    ierror = MPI_SUCCESS
    call mpi_copy( data1, data2, nsend, ierror )
  end subroutine mpi_allgather_real32_r1

  subroutine mpi_allgather_real64_r1( data1, nsend, sendtype,data2, nrecv, recvtype, comm, ierror )
    implicit none
    integer nsend, nrecv, sendtype, recvtype, comm, ierror
    real(c_double) :: data1(nsend), data2(nsend)
    ierror = MPI_SUCCESS
    call mpi_copy( data1, data2, nsend, ierror )
  end subroutine mpi_allgather_real64_r1

  !*****************************************************************************80
  !
  !! MPI_ALLGATHERV gathers data from all the processes in a communicator.
  !
  !*****************************************************************************80
  
  subroutine mpi_allgatherv_int_r1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
    implicit none
    integer :: nsend
    integer :: data1(nsend), data2(nsend)
    integer :: comm, ierror, ndispls, nrecv, recvtype, sendtype
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_allgatherv_int_r1
  
  subroutine mpi_allgatherv_real32_r1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
    implicit none
    integer :: nsend
    real :: data1(nsend), data2(nsend)
    integer :: comm, ierror, ndispls, nrecv, recvtype, sendtype
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_allgatherv_real32_r1
  
  subroutine mpi_allgatherv_real64_r1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, comm, ierror )
    implicit none
    integer :: nsend
    real(c_double) :: data1(nsend), data2(nsend)
    integer :: comm, ierror, ndispls, nrecv, recvtype, sendtype
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_allgatherv_real64_r1
  
  !*****************************************************************************80
  !
  !! MPI_GATHERV gathers data from all the processes in a communicator.
  !
  !*****************************************************************************80
  
  subroutine mpi_gatherv_int_r1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, root, comm, ierror )
    implicit none
    integer :: nsend
    integer :: data1(nsend), data2(nsend)
    integer :: comm, ierror, ndispls(:), nrecv(:), recvtype, sendtype, root
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_gatherv_int_r1
  
  subroutine mpi_gatherv_real32_r1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, root, comm, ierror )
    implicit none    
    integer :: nsend
    real :: data1(nsend), data2(nsend)
    integer :: comm, ierror, ndispls(:), nrecv(:), recvtype, sendtype, root
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_gatherv_real32_r1
  
  subroutine mpi_gatherv_real64_r1 ( data1, nsend, sendtype, data2, nrecv, ndispls, recvtype, root, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: nsend
    real(c_double) :: data1(nsend), data2(nsend)
    integer :: comm, ierror, ndispls(:), nrecv(:), recvtype, sendtype, root
    
    ierror = MPI_SUCCESS
    
    call mpi_copy( data1, data2, nsend, ierror )
    
  end subroutine mpi_gatherv_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_ALLTOALL carries out a reduction operation.
  !
  !*****************************************************************************
  
  subroutine mpi_alltoall_int_r1 ( data1, nsend, sendtype, data2, nrecv, recvtype, comm, ierror )
    integer :: nsend, nrecv
    integer :: data1(nsend), data2(nrecv)
    integer :: comm, sendtype, recvtype, ierror
    data2 = data1
  end subroutine mpi_alltoall_int_r1
  
  subroutine mpi_alltoall_real32_r1 ( data1, nsend, sendtype, data2, nrecv, recvtype, comm, ierror )
    integer :: nsend, nrecv
    real(c_float) :: data1(nsend), data2(nrecv)
    integer :: comm, sendtype, recvtype, ierror
    data2(:) = data1(:)
  end subroutine mpi_alltoall_real32_r1

  subroutine mpi_alltoall_real64_r1 ( data1, nsend, sendtype, data2, nrecv, recvtype, comm, ierror )
    integer :: nsend, nrecv
    real(c_double) :: data1(nsend), data2(nrecv)
    integer :: comm, sendtype, recvtype, ierror
    data2(:) = data1(:)
  end subroutine mpi_alltoall_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_ALLTOALL carries out a reduction operation.
  !
  !*****************************************************************************
  
  subroutine mpi_alltoallv_int_r1 ( data1, nsend, senddispls, sendtype, data2, nrecv, recvdispls, recvtype, comm, ierror )
    integer :: nsend(:), nrecv(:), senddispls(:), recvdispls(:)
    integer :: data1(:), data2(:)
    integer :: comm, sendtype, recvtype, ierror
    data2(:) = data1(:)
  end subroutine mpi_alltoallv_int_r1
  
  subroutine mpi_alltoallv_real32_r1 ( data1, nsend, senddispls, sendtype, data2, nrecv, recvdispls, recvtype, comm, ierror )
    integer :: nsend(:), nrecv(:), senddispls(:), recvdispls(:)
    real(c_float) :: data1(:), data2(:)
    integer :: comm, sendtype, recvtype, ierror
    data2(:) = data1(:)
  end subroutine mpi_alltoallv_real32_r1
    
  subroutine mpi_alltoallv_real64_r1 ( data1, nsend, senddispls, sendtype, data2, nrecv, recvdispls, recvtype, comm, ierror )
    integer :: nsend(:), nrecv(:), senddispls(:), recvdispls(:)
    real(c_double) :: data1(:), data2(:)
    integer :: comm, sendtype, recvtype, ierror
    data2(:) = data1(:)
  end subroutine mpi_alltoallv_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_ALLREDUCE carries out a reduction operation.
  !
  !*****************************************************************************
  
  subroutine mpi_allreduce_int_r0( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1, data2
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
    if( data1 .ne. MPI_IN_PLACE ) data2 = data1    
  end subroutine mpi_allreduce_int_r0
  
  subroutine mpi_allreduce_int_r1( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1(n), data2(n)
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
    data2 = data1    
  end subroutine mpi_allreduce_int_r1
  
  subroutine mpi_allreduce_int_r1_inplace( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1, data2(n)
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
  end subroutine mpi_allreduce_int_r1_inplace
  
  subroutine mpi_allreduce_real32_r0( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    real(c_float) :: data1, data2
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
    data2 = data1
  end subroutine mpi_allreduce_real32_r0

  subroutine mpi_allreduce_real64_r0( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    real(c_double) :: data1, data2
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
    data2 = data1
  end subroutine mpi_allreduce_real64_r0

  subroutine mpi_allreduce_real32_r0_inplace( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1
    real(c_float) :: data2
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
  end subroutine mpi_allreduce_real32_r0_inplace

  subroutine mpi_allreduce_real64_r0_inplace( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1
    real(c_double) :: data2
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
  end subroutine mpi_allreduce_real64_r0_inplace

  
  subroutine mpi_allreduce_real32_r1( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    real(c_float):: data1(n), data2(n)
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
    data2 = data1
  end subroutine mpi_allreduce_real32_r1

  subroutine mpi_allreduce_real64_r1( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    real(c_double)    :: data1(n), data2(n)
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
    data2 = data1
  end subroutine mpi_allreduce_real64_r1
  
  subroutine mpi_allreduce_real32_r1_inplace( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1
    real(c_float)    :: data2(n)
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS    
  end subroutine mpi_allreduce_real32_r1_inplace

  subroutine mpi_allreduce_real64_r1_inplace( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1
    real(c_double)    :: data2(n)
    integer :: comm, datatype, ierror, operation
    ierror = MPI_SUCCESS
  end subroutine mpi_allreduce_real64_r1_inplace

  !*****************************************************************************
  !
  ! MPI_BARRIER forces processes within a communicator to wait together.
  !
  !*****************************************************************************
  
  subroutine mpi_barrier ( comm, ierror )
    implicit none
    
    integer comm, ierror
    integer, parameter :: MPI_FAILURE = 1
    integer, parameter :: MPI_SUCCESS = 0
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_barrier
  

  !*****************************************************************************
  !
  ! MPI_REDUCE
  !
  !*****************************************************************************
  
  subroutine mpi_reduce_int_r0( data1, data2, n, datatype, operation, receiver, comm, ierror )
    implicit none
    integer :: n
    integer :: data1, data2
    integer :: datatype, operation, receiver, comm, ierror
    
    ierror = MPI_SUCCESS
    
    data2 = data1
    
    return
  end subroutine mpi_reduce_int_r0
  
  subroutine mpi_reduce_int_r1( data1, data2, n, datatype, operation, receiver, comm, ierror )
    implicit none
    integer :: n
    integer :: data1(n), data2(n)
    integer :: datatype, operation, receiver, comm, ierror
    
    ierror = MPI_SUCCESS

    data2 = data1
    
    return
  end subroutine mpi_reduce_int_r1
  
  subroutine mpi_reduce_real32_r0( data1, data2, n, datatype, operation, receiver, comm, ierror )
    implicit none
    integer :: n
    real(c_float) :: data1, data2
    integer :: datatype, operation, receiver, comm, ierror
    ierror = MPI_SUCCESS
    data2 = data1
    return
  end subroutine mpi_reduce_real32_r0
  
  subroutine mpi_reduce_real64_r0( data1, data2, n, datatype, operation, receiver, comm, ierror )
    implicit none
    integer :: n
    real(c_double) :: data1, data2
    integer :: datatype, operation, receiver, comm, ierror
    ierror = MPI_SUCCESS
    data2 = data1
    return
  end subroutine mpi_reduce_real64_r0

  subroutine mpi_reduce_real32_r1 ( data1, data2, n, datatype, operation, receiver, comm, ierror )
    implicit none
    integer :: n
    real(c_float) :: data1(n), data2(n)
    integer :: datatype, operation, receiver, comm, ierror
    ierror = MPI_SUCCESS
    data2 = data1
    return
  end subroutine mpi_reduce_real32_r1
  
  subroutine mpi_reduce_real64_r1 ( data1, data2, n, datatype, operation, receiver, comm, ierror )
    implicit none
    integer :: n
    real(c_double) :: data1(n), data2(n)
    integer :: datatype, operation, receiver, comm, ierror
    ierror = MPI_SUCCESS
    data2 = data1
    return
  end subroutine mpi_reduce_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_REDUCE_SCATTER collects a message of the same length from each process.
  !
  !*****************************************************************************
  
  subroutine mpi_reduce_scatter_int_r1( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    integer :: data1(n), data2(n)
    integer comm, datatype, ierror, operation
    
    ierror = MPI_SUCCESS
    
    data2(1:n) = data1(1:n)
    
  end subroutine mpi_reduce_scatter_int_r1
  
  subroutine mpi_reduce_scatter_real32_r1( data1, data2, n, datatype, operation, comm, ierror )
    implicit none
    integer :: n
    real :: data1(n), data2(n)
    integer comm, datatype, ierror, operation
    
    ierror = MPI_SUCCESS
    
    data2(1:n) = data1(1:n)
    
  end subroutine mpi_reduce_scatter_real32_r1
  
  
  
  !*****************************************************************************
  !
  ! MPI_WAITALL waits until all I/O requests have completed.
  !
  !*****************************************************************************
  
  subroutine mpi_waitall ( icount, irequest, istatus, ierror )
    implicit none
    
    integer icount, ierror
    integer irequest(icount), istatus(icount)
    integer, parameter :: MPI_FAILURE = 1
    integer, parameter :: MPI_SUCCESS = 0
    
    ierror = MPI_SUCCESS
  end subroutine mpi_waitall
  
  !*****************************************************************************
  !
  ! MPI_ISEND sends data from one process to another using nonblocking transmission.
  !
  !*****************************************************************************

  subroutine mpi_isend_int_r0 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    integer :: data
    integer :: comm, datatype, ierror, iproc, itag, request
    
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_isend_int_r0
  
  subroutine mpi_isend_int_r1 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    integer :: data(n)
    integer :: comm, datatype, ierror, iproc, itag, request
    
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_isend_int_r1
  
  subroutine mpi_isend_real32_r0 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    real :: data
    integer :: comm, datatype, ierror, iproc, itag, request
    
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_isend_real32_r0
  
  subroutine mpi_isend_real32_r1 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    real :: data(n)
    integer :: comm, datatype, ierror, iproc, itag, request
    
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_isend_real32_r1
  
  subroutine mpi_isend_real64_r0 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    real(c_double) :: data
    integer :: comm, datatype, ierror, iproc, itag, request
    
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_isend_real64_r0
  
  subroutine mpi_isend_real64_r1 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    real(c_double) :: data(n)
    integer :: comm, datatype, ierror, iproc, itag, request
    
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_isend_real64_r1


  !*****************************************************************************
  !
  ! MPI_IRECV sends data from one process to another using nonblocking transmission.
  !
  !*****************************************************************************

  subroutine mpi_irecv_int_r0 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    integer :: data
    integer :: comm, datatype, ierror, iproc, itag, request
    
    data = 0
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_irecv_int_r0
  
  subroutine mpi_irecv_int_r1 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    integer :: data(n)
    integer :: comm, datatype, ierror, iproc, itag, request
    
    data = 0
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_irecv_int_r1
  
  subroutine mpi_irecv_real32_r0 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    real :: data
    integer :: comm, datatype, ierror, iproc, itag, request
    
    data = 0
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_irecv_real32_r0
  
  subroutine mpi_irecv_real32_r1 ( data, n, datatype, iproc, itag, comm, request, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    real :: data(n)
    integer :: comm, datatype, ierror, iproc, itag, request
    data = 0
    request = 0
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_irecv_real32_r1
  
  subroutine mpi_irecv_real64_r0 ( data, n, datatype, iproc, itag, comm, request, ierror )
      implicit none
    
      !include "mpi_stubs_f90.h"
    
      integer :: n
      real(c_double) :: data
      integer :: comm, datatype, ierror, iproc, itag, request
    
      data = 0
      request = 0
      ierror = MPI_SUCCESS
    
      return
    end subroutine mpi_irecv_real64_r0
  
    subroutine mpi_irecv_real64_r1 ( data, n, datatype, iproc, itag, comm, request, ierror )
      implicit none
    
      !include "mpi_stubs_f90.h"
    
      integer :: n
      real(c_double) :: data(n)
      integer :: comm, datatype, ierror, iproc, itag, request
      data = 0
      request = 0
      ierror = MPI_SUCCESS
    
      return
    end subroutine mpi_irecv_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_SEND
  !
  !*****************************************************************************

  subroutine mpi_send_int_r0 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    integer :: buf
    integer :: datatype, dest, tag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_int_r0
  
  subroutine mpi_send_int_r1 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    integer :: buf(count)
    integer :: datatype, dest, tag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_int_r1
  
  subroutine mpi_send_int_r2 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    integer :: buf(:,:)
    integer :: datatype, dest, tag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_int_r2
  
  subroutine mpi_send_int_r2_fixedstatus ( buf, count, datatype, dest, tag, comm, status, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    integer :: buf(:,:)
    integer :: status(MPI_STATUS_SIZE)
    integer :: datatype, dest, tag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_int_r2_fixedstatus
  
  subroutine mpi_send_real32_r0 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    real :: buf
    integer :: datatype, dest, tag, comm, ierror

    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_real32_r0
  
  subroutine mpi_send_real32_r1 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    real :: buf(count)
    integer :: datatype, dest, tag, comm, ierror

    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_real32_r1
  
  subroutine mpi_send_real64_r0 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    real(c_double) :: buf
    integer :: datatype, dest, tag, comm, ierror

    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_real64_r0
  
  subroutine mpi_send_real64_r1 ( buf, count, datatype, dest, tag, comm, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: count
    real(c_double) :: buf(count)
    integer :: datatype, dest, tag, comm, ierror

    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_send_real64_r1
  
  !*****************************************************************************
  !
  ! MPI_RECV receives data from another process within a communicator.
  !
  !*****************************************************************************
  
  subroutine mpi_recv_int_r0_fixedstatus( data, n, datatype, iproc, itag, comm, istatus, ierror )
    implicit none
    
    !include "mpi_stubs_f90.h"
    
    integer :: n
    integer :: data, istatus(MPI_STATUS_SIZE)
    integer datatype, iproc, itag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_recv_int_r0_fixedstatus
  
  subroutine mpi_recv_int_r0( data, n, datatype, iproc, itag, comm, istatus, ierror )
    implicit none
    integer :: n
    integer :: data, istatus
    integer datatype, iproc, itag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_recv_int_r0
  
  subroutine mpi_recv_int_r1( data, n, datatype, iproc, itag, comm, istatus, ierror )
    implicit none
    integer :: n
    integer :: data(:), istatus(:)
    integer datatype, iproc, itag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_recv_int_r1
  
  subroutine mpi_recv_int_r2( data, n, datatype, iproc, itag, comm, istatus, ierror )
    implicit none
    integer :: n
    integer :: data(:,:), istatus(*)
    integer datatype, iproc, itag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_recv_int_r2

  subroutine mpi_recv_real32_r0( data, n, datatype, iproc, itag, comm, istatus, ierror )
    implicit none
    integer :: n
    real :: data
    integer :: istatus
    integer datatype, iproc, itag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_recv_real32_r0

  subroutine mpi_recv_real32_r1( data, n, datatype, iproc, itag, comm, istatus, ierror )
    implicit none
    integer :: n
    real :: data(n)
    integer :: istatus(n)
    integer datatype, iproc, itag, comm, ierror
    
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_recv_real32_r1

  !*****************************************************************************
  !
  ! MPI_PACK
  !
  !*****************************************************************************
  
  subroutine mpi_pack_int_r0 ( inbuf, incount, datatype, outbuf, outcount, position, comm, ierror )
    
    implicit none
    
    integer incount, outcount
    integer inbuf
    integer outbuf(outcount)
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_pack_int_r0

  subroutine mpi_pack_int_r1 ( inbuf, incount, datatype, outbuf, outcount, position, comm, ierror )
    implicit none
    integer incount, outcount
    integer inbuf(incount)
    integer outbuf(outcount)
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    
    return
  end subroutine mpi_pack_int_r1

  subroutine mpi_pack_real32_r0 ( inbuf, incount, datatype, outbuf, outcount, position, comm, ierror )
    implicit none
    integer incount, outcount
    real inbuf
    integer outbuf(outcount)
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    return
  end subroutine mpi_pack_real32_r0

  subroutine mpi_pack_real32_r1 ( inbuf, incount, datatype, outbuf, outcount, position, comm, ierror )
    implicit none
    integer incount, outcount
    real inbuf(incount)
    integer outbuf(outcount)
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    return
  end subroutine mpi_pack_real32_r1

  !*****************************************************************************
  !
  ! MPI_UNPACK
  !
  !*****************************************************************************
  
  subroutine mpi_unpack_int_r0( inbuf, insize, position, outbuf, outcount, datatype, comm, ierror )
    
    implicit none
    
    integer insize, outcount
    integer inbuf(insize)
    integer outbuf
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    
    return
    
  end subroutine mpi_unpack_int_r0

  subroutine mpi_unpack_int_r1( inbuf, insize, position, outbuf, outcount, datatype, comm, ierror )
    
    implicit none
    
    integer insize, outcount
    integer inbuf(insize)
    integer outbuf(outcount)
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    
    return
    
  end subroutine mpi_unpack_int_r1

  subroutine mpi_unpack_real32_r0( inbuf, insize, position, outbuf, outcount, datatype, comm, ierror )
    
    implicit none
    
    integer insize, outcount
    integer inbuf(insize)
    real outbuf
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    return
    
  end subroutine mpi_unpack_real32_r0

  subroutine mpi_unpack_real32_r1( inbuf, insize, position, outbuf, outcount, datatype, comm, ierror )
    implicit none
    integer insize, outcount
    integer inbuf(insize)
    real outbuf(outcount)
    integer datatype
    integer position
    integer comm
    integer ierror
    ierror = MPI_SUCCESS
    return
  end subroutine mpi_unpack_real32_r1

!*****************************************************************************80
!
!! MPI_BCAST broadcasts data from one process to all others.
!
!*****************************************************************************80

subroutine mpi_bcast_int_r0( data, n, datatype, node, comm, ierror )
  implicit none

  integer n

  integer comm
  integer data
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node

  ierror = MPI_SUCCESS

  return
end subroutine mpi_bcast_int_r0

subroutine mpi_bcast_int_r1( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node

  ierror = MPI_SUCCESS

  return
end subroutine mpi_bcast_int_r1

subroutine mpi_bcast_int_r2( data, n, datatype, node, comm, ierror )
  implicit none

  integer n

  integer comm
  integer data(:,:)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node

  ierror = MPI_SUCCESS

  return
end subroutine mpi_bcast_int_r2

subroutine mpi_bcast_real32_r1( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  real(c_float) data(n)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node

  ierror = MPI_SUCCESS

  return
end subroutine mpi_bcast_real32_r1

subroutine mpi_bcast_real32_r0( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  real(c_float) data
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node
  ierror = MPI_SUCCESS
  return
end subroutine mpi_bcast_real32_r0

subroutine mpi_bcast_real64_r0( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  real(c_double) data
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node
  ierror = MPI_SUCCESS
  return
end subroutine mpi_bcast_real64_r0

subroutine mpi_bcast_real64_r1( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  real(c_double) data(n)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node
  ierror = MPI_SUCCESS
  return
end subroutine mpi_bcast_real64_r1

subroutine mpi_bcast_real32_r2( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  real(c_float) data(:,:)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node
  ierror = MPI_SUCCESS
  return
end subroutine mpi_bcast_real32_r2

subroutine mpi_bcast_real64_r2( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  real(c_double) data(:,:)
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node
  ierror = MPI_SUCCESS
  return
end subroutine mpi_bcast_real64_r2

subroutine mpi_bcast_char( data, n, datatype, node, comm, ierror )
  implicit none
  integer n
  integer comm
  character(LEN=n) data
  integer datatype
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer node
  ierror = MPI_SUCCESS
  return
end subroutine mpi_bcast_char


subroutine mpi_probe(data,itag,comm,istatus,ierror)
  implicit none
  integer comm
  integer data
  integer itag
  integer istatus(*)
  integer ierror
  ierror = MPI_SUCCESS
  return
end subroutine mpi_probe

subroutine mpi_iprobe(data,itag,comm,flag,istatus,ierror)
  implicit none
  integer comm
  integer data
  integer itag
  integer istatus
  integer ierror
  logical flag
  ierror = MPI_SUCCESS
  flag   = .false.
  return
end subroutine mpi_iprobe


subroutine mpi_test(request,flag,istatus,ierror)
  implicit none
  integer request
  logical flag
  integer istatus
  integer ierror
  ierror = MPI_SUCCESS
  return
end subroutine mpi_test















!*****************************************************************************80
!
!! MPI_BSEND sends data from one process to another, using buffering.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Modified:
!
!    06 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, datatype DATA(N), the data to be sent.
!
!    Input, integer N, the number of data items to send.
!
!    Input, integer DATAYTPE, the MPI code for the datatype.
!
!    Input, integer IPROC, the rank of the process within the communicator
!    that is to receive the message.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_bsend ( data, n, datatype, iproc, itag, comm, ierror )
  implicit none
  integer n
  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_BSEND - Error!'
  write ( *, '(a)' )  '  Should not send message to self.'

  return
end subroutine mpi_bsend

!*****************************************************************************80
!
!! MPI_CART_CREATE creates a communicator for a Cartesian topology.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_cart_create ( comm, ndims, dims, periods, reorder, comm_cart, &
  ierror )

  implicit none

  integer ndims

  integer comm
  integer comm_cart
  integer dims(*)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  logical periods(*)
  logical reorder

  ierror = MPI_SUCCESS

  return
end subroutine mpi_cart_create

!*****************************************************************************80
!
!! MPI_CART_GET returns the "Cartesian coordinates" of the calling process.
!
!  Discussion:
!
!    Set all coordinates to 0.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_cart_get ( comm, ndims, dims, periods, coords, ierror )
  implicit none

  integer ndims

  integer comm
  integer coords(*)
  integer dims(*)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  logical periods(*)

  ierror = MPI_SUCCESS

  coords(1:ndims) = 0

  return
end subroutine mpi_cart_get


!*****************************************************************************80
!
!! MPI_CART_SHIFT finds the destination and source for Cartesian shifts.
!
!  Discussion:
!
!    Set ISOURCE = IDEST = SELF = 0.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_cart_shift ( comm, idir, idisp, isource, idest, ierror )
  implicit none

  integer comm
  integer idest
  integer idir
  integer idisp
  integer ierror
  integer isource
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS
  isource = 0
  idest = 0

  return
end subroutine mpi_cart_shift

!*****************************************************************************80
!
!! MPI_COMM_DUP duplicates a communicator.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_comm_dup ( comm, comm_out, ierror )
  implicit none

  integer comm
  integer comm_out
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine mpi_comm_dup


!*****************************************************************************80
!
!! MPI_COMM_FREE "frees" a communicator.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_comm_free ( comm, ierror )
  implicit none

  integer comm
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine mpi_comm_free


!*****************************************************************************80
!
!! MPI_COMM_RANK reports the rank of the calling process.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_comm_rank ( comm, me, ierror )
  implicit none

  integer comm
  integer ierror
  integer me
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS
  me = 0

  return
end subroutine mpi_comm_rank


!*****************************************************************************80
!
!! MPI_COMM_SIZE reports the number of processes in a communicator.
!
!  Discussion:
!
!    The routine simply returns NPROCS = 1.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_comm_size ( comm, nprocs, ierror )
  implicit none

  integer comm
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0
  integer nprocs

  ierror = MPI_SUCCESS
  nprocs = 1

  return
end subroutine mpi_comm_size


!*****************************************************************************80
!
!! MPI_COMM_SPLIT splits up a communicator based on a key.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_comm_split ( comm, icolor, ikey, comm_new, ierror )
  implicit none

  integer comm
  integer comm_new
  integer icolor
  integer ierror
  integer ikey
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine mpi_comm_split


!*****************************************************************************80
!
!! MPI_COPY_DOUBLE copies a double precision vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, double precision DATA1(N), the data to be copied.
!
!    Output, double precision DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_copy_real64_r1 ( data1, data2, n, ierror )
  implicit none

  integer n

  real(c_double) data1(n)
  real(c_double) data2(n)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine mpi_copy_real64_r1


!*****************************************************************************80
!
!! MPI_COPY_INTEGER copies an integer vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer DATA1(N), the data to be copied.
!
!    Output, integer DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_copy_int_r1 ( data1, data2, n, ierror )
  implicit none

  integer n

  integer data1(n)
  integer data2(n)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine mpi_copy_int_r1


!*****************************************************************************80
!
!! MPI_COPY_REAL copies a real vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, real DATA1(N), the data to be copied.
!
!    Output, real DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_copy_real32_r1 ( data1, data2, n, ierror )
  implicit none

  integer n

  real data1(n)
  real data2(n)
  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine mpi_copy_real32_r1


!*****************************************************************************80
!
!! MPI_FINALIZE shuts down the MPI library.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_finalize ( ierror )
  implicit none

  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine mpi_finalize

subroutine mpi_finalized ( flag, ierror )
  implicit none

  logical ::flag
  integer :: ierror

  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  flag = .True.
  ierror = MPI_SUCCESS

  return
end subroutine mpi_finalized

!*****************************************************************************80
!
!! MPI_GET_COUNT reports the actual number of items transmitted.
!
!  Discussion:
!
!    Warn against querying message from self, since no data copy is done.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_get_count ( istatus, datatype, icount, ierror )
  implicit none

  integer datatype
  integer icount
  integer ierror
  integer istatus(icount)
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_GET_COUNT - Error!'
  write ( *, '(a)' ) '  Should not query message from self.'

  return
end subroutine mpi_get_count


!*****************************************************************************80
!
!! MPI_INIT initializes the MPI library.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_init ( ierror )
  implicit none

  integer ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_SUCCESS

  return
end subroutine mpi_init

subroutine mpi_initialized ( flag, ierror )
  implicit none

  logical ::flag
  integer :: ierror
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  flag = .True.
  ierror = MPI_SUCCESS

  return
end subroutine mpi_initialized

!*****************************************************************************80
!
!! MPI_IRECV receives data from another process.
!
!  Discussion:
!
!    Warn against receiving message from self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_irecv_gen ( data, n, datatype, iproc, itag, comm, irequest, ierror )
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer irequest
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_IRECV - Error!'
  write ( *, '(a)' ) '  Should not recv message from self.'

  return
end subroutine mpi_irecv_gen


!*****************************************************************************80
!
!! MPI_RSEND "ready sends" data from one process to another.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type, 
!    but using the other types should generally not cause a problem.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_rsend ( data, n, datatype, iproc, itag, comm, ierror )
  implicit none

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_RSEND - Error!'
  write ( *, '(a)' ) '  Should not send message to self.'

  return
end subroutine mpi_rsend

!*****************************************************************************80
!
!! MPI_WAIT waits for an I/O request to complete.
!
!  Discussion:
!
!    Warn against waiting on message from self, since no data copy is done.
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_wait ( irequest, istatus, ierror )
  implicit none

  integer ierror
  integer irequest
  integer istatus(MPI_STATUS_SIZE)
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_WAIT - Error!'
  write ( *, '(a)' ) '  Should not wait on message from self.'

  return
end subroutine mpi_wait


!*****************************************************************************80
!
!! MPI_WAITANY waits until one I/O requests has completed.
!
!  Discussion:
!
!    Warn against waiting on message from self, since no data copy is done.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
subroutine mpi_waitany ( icount, array_of_requests, index, istatus, ierror )
  implicit none

  integer array_of_requests(*)
  integer icount
  integer ierror
  integer index
  integer istatus
  integer, parameter :: MPI_FAILURE = 1
  integer, parameter :: MPI_SUCCESS = 0

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_WAITANY - Error!'
  write ( *, '(a)' ) '  Should not wait on message from self.'

  return
end subroutine mpi_waitany

!*****************************************************************************80
!
!! MPI_WTICK returns the number of seconds per clock tick.
!
!  Discussion:
!
!    The value returned here is simply a dummy value.
!
!  Modified:
!
!    04 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MPI_WTICK, the number of seconds per clock tick.
!
function mpi_wtick ( )
  implicit none

  real ( kind = 8 ) mpi_wtick
  
  mpi_wtick = 1.0D+00
  
  return
end function mpi_wtick


!*****************************************************************************80
!
!! MPI_WTIME returns the elapsed wall clock time.
!
!  Modified:
!
!    26 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MPI_WTIME, the elapsed wall clock time.
!
function mpi_wtime ( )
  implicit none

  integer count
  integer count_max
  integer count_rate
  real ( kind = 8 ) mpi_wtime

  call system_clock ( count, count_rate, count_max )
  
  mpi_wtime = real ( count, kind = 8 ) / real ( count_rate, kind = 8 )
  
  return
end function mpi_wtime

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
subroutine timestring ( string )
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestring

end module atlas_mpistubs_module
