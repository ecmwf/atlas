! =====================================================================
! parallel_module
! This module contains routines for MPI communication
! =====================================================================

module parallel_module
  use mpi
  implicit none
  private

public myproc, nproc, Comm_type, parallel_init, parallel_finalise

  integer :: ierr
  integer :: myproc = -1
  integer :: nproc = -1


  type :: Comm_type
  private
    integer, dimension(:), allocatable :: sendcounts
    integer, dimension(:), allocatable :: senddispl
    integer, dimension(:), allocatable :: recvcounts
    integer, dimension(:), allocatable :: recvdispl
    integer, dimension(:), allocatable :: sendmap
    integer, dimension(:), allocatable :: recvmap
    integer :: sendcnt, recvcnt
  contains
    procedure, public :: setup => setup_comm
    procedure, pass :: synchronise_real8_rank1
    procedure, pass :: synchronise_real8_rank2
    procedure, pass :: synchronise_integer_rank1
    procedure, pass :: synchronise_integer_rank2
    generic, public :: synchronise => &
      & synchronise_real8_rank1, &
      & synchronise_real8_rank2, &
      & synchronise_integer_rank1, &
      & synchronise_integer_rank2

  end type Comm_type


contains

  subroutine parallel_init()
    call MPI_INIT (ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD,myproc,ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
  end subroutine parallel_init

  subroutine parallel_finalise()
    call MPI_FINALIZE (ierr)
  end subroutine parallel_finalise

  subroutine setup_comm(comm,proc,glb_idx)
    class(Comm_type), intent(inout) :: comm
    integer, dimension(:), intent(in) :: glb_idx
    integer, dimension(:), intent(in) :: proc
    integer :: max_glb_idx = -1
    integer :: jnode, jproc, ierr, cnt(0:nproc-1), nb_nodes
    integer, dimension(:), allocatable :: send_requests, recv_requests
    integer, dimension(:), allocatable :: map_glb_to_loc

    allocate( comm%sendcounts(0:nproc-1) )
    allocate( comm%recvcounts(0:nproc-1) )
    allocate( comm%senddispl(0:nproc-1) )
    allocate( comm%recvdispl(0:nproc-1) )
  
    nb_nodes = size(glb_idx)

    do jnode=1,nb_nodes
      max_glb_idx = max( max_glb_idx, glb_idx(jnode) )
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, max_glb_idx, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    allocate( map_glb_to_loc(max_glb_idx) )

    comm%recvcounts(:) = 0
    map_glb_to_loc(:) = -1
    do jnode=1,nb_nodes
      if (proc(jnode) .ne. myproc) then
        comm%recvcounts(proc(jnode)) = comm%recvcounts(proc(jnode))+1
      end if
      map_glb_to_loc(glb_idx(jnode)) = jnode
    end do

    call MPI_ALLTOALL( comm%recvcounts, 1, MPI_INTEGER, &
                     & comm%sendcounts, 1, MPI_INTEGER, &
                     & MPI_COMM_WORLD,ierr )

    comm%recvcnt = sum(comm%recvcounts)
    comm%sendcnt = sum(comm%sendcounts)
    allocate( send_requests(comm%recvcnt) )
    allocate( recv_requests(comm%sendcnt) )
    allocate( comm%sendmap(comm%sendcnt) )
    allocate( comm%recvmap(comm%recvcnt) )

    comm%recvdispl(0)=0
    do jproc=1,nproc-1
      comm%recvdispl(jproc)=comm%recvcounts(jproc-1)+comm%recvdispl(jproc-1)
    end do

    comm%senddispl(0)=0
    do jproc=1,nproc-1
      comm%senddispl(jproc)=comm%sendcounts(jproc-1)+comm%senddispl(jproc-1)
    end do

    cnt(:)=0
    do jnode=1,nb_nodes
      if (proc(jnode) .ne. myproc ) then
        send_requests( comm%recvdispl(proc(jnode))+1+cnt(proc(jnode)) ) = glb_idx(jnode)
        comm%recvmap( comm%recvdispl(proc(jnode))+1+cnt(proc(jnode)) ) = jnode 
        cnt(proc(jnode)) = cnt(proc(jnode))+1
      endif
    end do

    call MPI_ALLTOALLV( send_requests,comm%recvcounts,comm%recvdispl,MPI_INTEGER, &
                      & recv_requests,comm%sendcounts,comm%senddispl,MPI_INTEGER, &
                      & MPI_COMM_WORLD,ierr )
  
    do jnode=1,comm%sendcnt
      comm%sendmap(jnode) = map_glb_to_loc( recv_requests(jnode) )
    end do

  end subroutine setup_comm



  subroutine synchronise_real8_rank1(comm,field)
    class(Comm_type), intent(inout) :: comm
    real*8, dimension(:), intent(inout) :: field
    real*8 :: sendbuffer(comm%sendcnt)
    real*8 :: recvbuffer(comm%recvcnt)
    integer :: jnode

    ! Pack
    do jnode=1,comm%sendcnt
      sendbuffer(jnode) = field(comm%sendmap(jnode))
    end do

    ! Communicate
    call MPI_ALLTOALLV( sendbuffer,comm%sendcounts,comm%senddispl,MPI_DOUBLE_PRECISION, &
                      & recvbuffer,comm%recvcounts,comm%recvdispl,MPI_DOUBLE_PRECISION, &
                      & MPI_COMM_WORLD,ierr)
  
    ! Unpack
    do jnode=1,comm%recvcnt
      field( comm%recvmap(jnode) ) = recvbuffer(jnode)
    end do

  end subroutine synchronise_real8_rank1


  subroutine synchronise_real8_rank2(comm,field)
    class(Comm_type), intent(inout) :: comm
    real*8, dimension(:,:), intent(inout) :: field
    real*8 :: sendbuffer(comm%sendcnt)
    real*8 :: recvbuffer(comm%recvcnt)
    integer :: jnode
    integer :: jcol, ncols

    ncols = size(field,2)
    do jcol=1,ncols

      ! Pack
      do jnode=1,comm%sendcnt
        sendbuffer(jnode) = field(comm%sendmap(jnode),jcol)
      end do

      ! Communicate
      call MPI_ALLTOALLV( sendbuffer,comm%sendcounts,comm%senddispl,MPI_DOUBLE_PRECISION, &
                        & recvbuffer,comm%recvcounts,comm%recvdispl,MPI_DOUBLE_PRECISION, &
                        & MPI_COMM_WORLD,ierr)
    
      ! Unpack
      do jnode=1,comm%recvcnt
        field( comm%recvmap(jnode), jcol ) = recvbuffer(jnode)
      end do
    end do
  end subroutine synchronise_real8_rank2


  subroutine synchronise_integer_rank1(comm,field)
    class(Comm_type), intent(inout) :: comm
    integer, dimension(:), intent(inout) :: field
    integer :: sendbuffer(comm%sendcnt)
    integer :: recvbuffer(comm%recvcnt)
    integer :: jnode

    ! Pack
    do jnode=1,comm%sendcnt
      sendbuffer(jnode) = field(comm%sendmap(jnode))
    end do

    ! Communicate
    call MPI_ALLTOALLV( sendbuffer,comm%sendcounts,comm%senddispl,MPI_INTEGER, &
                      & recvbuffer,comm%recvcounts,comm%recvdispl,MPI_INTEGER, &
                      & MPI_COMM_WORLD,ierr)
  
    ! Unpack
    do jnode=1,comm%recvcnt
      field( comm%recvmap(jnode) ) = recvbuffer(jnode)
    end do
  end subroutine synchronise_integer_rank1


  subroutine synchronise_integer_rank2(comm,field)
    class(Comm_type), intent(inout) :: comm
    integer, dimension(:,:), intent(inout) :: field
    integer :: sendbuffer(comm%sendcnt)
    integer :: recvbuffer(comm%recvcnt)
    integer :: jnode
    integer :: jcol, ncols

    ncols = size(field,2)
    do jcol=1,ncols

      ! Pack
      do jnode=1,comm%sendcnt
        sendbuffer(jnode) = field(comm%sendmap(jnode),jcol)
      end do

      ! Communicate
      call MPI_ALLTOALLV( sendbuffer,comm%sendcounts,comm%senddispl,MPI_INTEGER, &
                        & recvbuffer,comm%recvcounts,comm%recvdispl,MPI_INTEGER, &
                        & MPI_COMM_WORLD,ierr)
    
      ! Unpack
      do jnode=1,comm%recvcnt
        field( comm%recvmap(jnode), jcol ) = recvbuffer(jnode)
      end do
    end do
  end subroutine synchronise_integer_rank2

end module parallel_module
