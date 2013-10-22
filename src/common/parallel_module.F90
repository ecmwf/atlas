! =====================================================================
! parallel_module
! This module contains routines for MPI communication
! =====================================================================

module parallel_module
  use mpi
  implicit none
  private

public myproc, nproc
public parallel_init, parallel_finalise, parallel_barrier
public Comm_type

  integer :: ierr
  integer :: myproc = -1
  integer :: nproc = -1


  type :: Comm_type
  private
    integer                            :: sync_sendcnt
    integer                            :: sync_recvcnt
    integer, dimension(:), allocatable :: sync_sendcounts
    integer, dimension(:), allocatable :: sync_senddispls
    integer, dimension(:), allocatable :: sync_recvcounts
    integer, dimension(:), allocatable :: sync_recvdispls
    integer, dimension(:), allocatable :: sync_sendmap
    integer, dimension(:), allocatable :: sync_recvmap

    integer, public :: gather_root = 0
    integer, public :: gather_sendcnt
    integer, public :: gather_recvcnt
    integer, public, dimension(:), allocatable :: gather_recvcounts
    integer, public, dimension(:), allocatable :: gather_recvdispls
    integer, public, dimension(:), allocatable :: gather_sendmap
    integer, public, dimension(:), allocatable :: gather_recvmap

  contains
    procedure, pass :: setup_gather  => setup_gather
    procedure, pass :: setup_sync    => setup_sync
    procedure, pass, public :: setup => setup_comm

    procedure, pass :: synchronise_real8_rank1
    procedure, pass :: synchronise_real8_rank2
    procedure, pass :: synchronise_integer_rank1
    procedure, pass :: synchronise_integer_rank2
    generic, public :: synchronise => &
      & synchronise_real8_rank1, &
      & synchronise_real8_rank2, &
      & synchronise_integer_rank1, &
      & synchronise_integer_rank2

    procedure, pass :: glb_field_size
    procedure, pass :: gather_real8_rank1
    procedure, pass :: gather_real8_rank2
    generic, public :: gather => &
      & gather_real8_rank1, &
      & gather_real8_rank2

  end type Comm_type


contains

  subroutine parallel_init()
    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myproc, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc,  ierr )
  end subroutine parallel_init

  subroutine parallel_finalise()
    call MPI_FINALIZE (ierr)
  end subroutine parallel_finalise

  subroutine parallel_barrier()
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  end subroutine

  subroutine setup_sync(comm,proc,glb_idx)
    class(Comm_type), intent(inout) :: comm
    integer, dimension(:), intent(in) :: glb_idx
    integer, dimension(:), intent(in) :: proc
    integer :: max_glb_idx = -1
    integer :: jnode, jproc, ierr, cnt(0:nproc-1), nb_nodes
    integer, dimension(:), allocatable :: send_requests, recv_requests
    integer, dimension(:), allocatable :: map_glb_to_loc

    allocate( comm%sync_sendcounts(0:nproc-1) )
    allocate( comm%sync_recvcounts(0:nproc-1) )
    allocate( comm%sync_senddispls(0:nproc-1) )
    allocate( comm%sync_recvdispls(0:nproc-1) )

    nb_nodes = size(glb_idx)

    do jnode=1,nb_nodes
      max_glb_idx = max( max_glb_idx, glb_idx(jnode) )
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, max_glb_idx, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    allocate( map_glb_to_loc(max_glb_idx) )

    comm%sync_recvcounts(:) = 0
    map_glb_to_loc(:) = -1
    do jnode=1,nb_nodes
      if (proc(jnode) .ne. myproc) then
        comm%sync_recvcounts(proc(jnode)) = comm%sync_recvcounts(proc(jnode))+1
      end if
      map_glb_to_loc(glb_idx(jnode)) = jnode
    end do

    call MPI_ALLTOALL( comm%sync_recvcounts, 1, MPI_INTEGER, &
                     & comm%sync_sendcounts, 1, MPI_INTEGER, &
                     & MPI_COMM_WORLD,ierr )

    comm%sync_recvcnt = sum(comm%sync_recvcounts)
    comm%sync_sendcnt = sum(comm%sync_sendcounts)
    allocate( send_requests(comm%sync_recvcnt) )
    allocate( recv_requests(comm%sync_sendcnt) )
    allocate( comm%sync_sendmap(comm%sync_sendcnt) )
    allocate( comm%sync_recvmap(comm%sync_recvcnt) )

    comm%sync_recvdispls(0)=0
    do jproc=1,nproc-1
      comm%sync_recvdispls(jproc)=comm%sync_recvcounts(jproc-1)+comm%sync_recvdispls(jproc-1)
    end do

    comm%sync_senddispls(0)=0
    do jproc=1,nproc-1
      comm%sync_senddispls(jproc)=comm%sync_sendcounts(jproc-1)+comm%sync_senddispls(jproc-1)
    end do

    cnt(:)=0
    do jnode=1,nb_nodes
      if (proc(jnode) .ne. myproc ) then
        send_requests( comm%sync_recvdispls(proc(jnode))+1+cnt(proc(jnode)) ) = glb_idx(jnode)
        comm%sync_recvmap( comm%sync_recvdispls(proc(jnode))+1+cnt(proc(jnode)) ) = jnode 
        cnt(proc(jnode)) = cnt(proc(jnode))+1
      end if
    end do

    call MPI_ALLTOALLV( send_requests,comm%sync_recvcounts,comm%sync_recvdispls,MPI_INTEGER, &
                      & recv_requests,comm%sync_sendcounts,comm%sync_senddispls,MPI_INTEGER, &
                      & MPI_COMM_WORLD,ierr )
  
    do jnode=1,comm%sync_sendcnt
      comm%sync_sendmap(jnode) = map_glb_to_loc( recv_requests(jnode) )
    end do

  end subroutine setup_sync

  subroutine setup_gather(comm,proc,glb_idx)
    class(Comm_type), intent(inout)   :: comm
    integer, dimension(:), intent(in) :: glb_idx
    integer, dimension(:), intent(in) :: proc
    integer :: jnode, jproc, ierr, nb_nodes, idx_send
    integer, dimension(:), allocatable :: gather_send_glb_idx

    allocate( comm%gather_recvcounts(0:nproc-1) )
    allocate( comm%gather_recvdispls(0:nproc-1) )

    nb_nodes = size(glb_idx)

    comm%gather_sendcnt = 0
    do jnode=1,nb_nodes
      if (proc(jnode) .eq. myproc) then
        comm%gather_sendcnt = comm%gather_sendcnt + 1
      end if
    end do

    comm%gather_recvcounts(:) = 0
    comm%gather_recvdispls(:) = 0
    call MPI_GATHER( comm%gather_sendcnt, 1, MPI_INTEGER, &
                   & comm%gather_recvcounts, 1, MPI_INTEGER, &
                   & comm%gather_root, MPI_COMM_WORLD, ierr );

    if (myproc .eq. comm%gather_root) then
      do jproc=1,nproc-1
        comm%gather_recvdispls(jproc) = comm%gather_recvdispls(jproc-1) + comm%gather_recvcounts(jproc-1)
      end do
    end if
    comm%gather_recvcnt = sum(comm%gather_recvcounts)
    allocate( comm%gather_sendmap(comm%gather_sendcnt) )
    allocate( comm%gather_recvmap(comm%gather_recvcnt) )
    allocate( gather_send_glb_idx(comm%gather_sendcnt) )

    idx_send = 0
    do jnode=1,nb_nodes
      if (proc(jnode) .eq. myproc ) then
        idx_send = idx_send + 1
        comm%gather_sendmap(idx_send) = jnode
        gather_send_glb_idx(idx_send) = glb_idx(jnode)
      endif
    end do

    ! This assumes that the global indices are one contiguous numbering across all procs
    call MPI_GATHERV( gather_send_glb_idx, comm%gather_sendcnt, MPI_INTEGER, &
                    & comm%gather_recvmap, comm%gather_recvcounts, comm%gather_recvdispls, MPI_INTEGER, &
                    & comm%gather_root, MPI_COMM_WORLD, ierr )

  end subroutine setup_gather

  subroutine setup_comm(comm,proc,glb_idx)
    class(Comm_type), intent(inout)   :: comm
    integer, dimension(:), intent(in) :: glb_idx
    integer, dimension(:), intent(in) :: proc
    call comm%setup_gather(proc,glb_idx)
    call comm%setup_sync(proc,glb_idx)
  end subroutine setup_comm

  subroutine synchronise_real8_rank1(comm,field)
    class(Comm_type), intent(inout) :: comm
    real*8, dimension(:), intent(inout) :: field
    real*8 :: sendbuffer(comm%sync_sendcnt)
    real*8 :: recvbuffer(comm%sync_recvcnt)
    integer :: jnode

    ! Pack
    do jnode=1,comm%sync_sendcnt
      sendbuffer(jnode) = field(comm%sync_sendmap(jnode))
    end do

    ! Communicate
    call MPI_ALLTOALLV( sendbuffer,comm%sync_sendcounts,comm%sync_senddispls,MPI_DOUBLE_PRECISION, &
                      & recvbuffer,comm%sync_recvcounts,comm%sync_recvdispls,MPI_DOUBLE_PRECISION, &
                      & MPI_COMM_WORLD,ierr)
  
    ! Unpack
    do jnode=1,comm%sync_recvcnt
      field( comm%sync_recvmap(jnode) ) = recvbuffer(jnode)
    end do

  end subroutine synchronise_real8_rank1


  subroutine synchronise_real8_rank2(comm,field)
    class(Comm_type), intent(inout) :: comm
    real*8, dimension(:,:), intent(inout) :: field
    real*8 :: sendbuffer(comm%sync_sendcnt)
    real*8 :: recvbuffer(comm%sync_recvcnt)
    integer :: jnode
    integer :: jcol, ncols

    ncols = size(field,2)
    do jcol=1,ncols

      ! Pack
      do jnode=1,comm%sync_sendcnt
        sendbuffer(jnode) = field(comm%sync_sendmap(jnode),jcol)
      end do

      ! Communicate
      call MPI_ALLTOALLV( sendbuffer,comm%sync_sendcounts,comm%sync_senddispls,MPI_DOUBLE_PRECISION, &
                        & recvbuffer,comm%sync_recvcounts,comm%sync_recvdispls,MPI_DOUBLE_PRECISION, &
                        & MPI_COMM_WORLD,ierr)
    
      ! Unpack
      do jnode=1,comm%sync_recvcnt
        field( comm%sync_recvmap(jnode), jcol ) = recvbuffer(jnode)
      end do
    end do
  end subroutine synchronise_real8_rank2


  subroutine synchronise_integer_rank1(comm,field)
    class(Comm_type), intent(inout) :: comm
    integer, dimension(:), intent(inout) :: field
    integer :: sendbuffer(comm%sync_sendcnt)
    integer :: recvbuffer(comm%sync_recvcnt)
    integer :: jnode

    ! Pack
    do jnode=1,comm%sync_sendcnt
      sendbuffer(jnode) = field(comm%sync_sendmap(jnode))
    end do

    ! Communicate
    call MPI_ALLTOALLV( sendbuffer,comm%sync_sendcounts,comm%sync_senddispls,MPI_INTEGER, &
                      & recvbuffer,comm%sync_recvcounts,comm%sync_recvdispls,MPI_INTEGER, &
                      & MPI_COMM_WORLD,ierr)
  
    ! Unpack
    do jnode=1,comm%sync_recvcnt
      field( comm%sync_recvmap(jnode) ) = recvbuffer(jnode)
    end do
  end subroutine synchronise_integer_rank1


  subroutine synchronise_integer_rank2(comm,field)
    class(Comm_type), intent(inout) :: comm
    integer, dimension(:,:), intent(inout) :: field
    integer :: sendbuffer(comm%sync_sendcnt)
    integer :: recvbuffer(comm%sync_recvcnt)
    integer :: jnode
    integer :: jcol, ncols

    ncols = size(field,2)
    do jcol=1,ncols

      ! Pack
      do jnode=1,comm%sync_sendcnt
        sendbuffer(jnode) = field(comm%sync_sendmap(jnode),jcol)
      end do

      ! Communicate
      call MPI_ALLTOALLV( sendbuffer,comm%sync_sendcounts,comm%sync_senddispls,MPI_INTEGER, &
                        & recvbuffer,comm%sync_recvcounts,comm%sync_recvdispls,MPI_INTEGER, &
                        & MPI_COMM_WORLD,ierr)
    
      ! Unpack
      do jnode=1,comm%sync_recvcnt
        field( comm%sync_recvmap(jnode), jcol ) = recvbuffer(jnode)
      end do
    end do
  end subroutine synchronise_integer_rank2


  subroutine gather_real8_rank1(comm,loc_field,glb_field)
    class(Comm_type), intent(inout) :: comm
    real*8, dimension(:), intent(in) :: loc_field
    real*8, dimension(:), allocatable, intent(inout) :: glb_field
    real*8 :: sendbuffer(comm%gather_sendcnt)
    real*8 :: recvbuffer(comm%gather_recvcnt)
    integer :: jnode

    if( .not. allocated( glb_field) ) then
      allocate( glb_field( comm%glb_field_size() ) )
    end if

    ! Pack
    do jnode=1,comm%gather_sendcnt
      sendbuffer(jnode) = loc_field(comm%gather_sendmap(jnode))
    end do

    ! Communicate
    call MPI_GATHERV( sendbuffer,comm%gather_sendcnt,MPI_DOUBLE_PRECISION, &
                    & recvbuffer,comm%gather_recvcounts,comm%gather_recvdispls,MPI_DOUBLE_PRECISION, &
                    & comm%gather_root,MPI_COMM_WORLD,ierr )
  
    ! Unpack
    do jnode=1,comm%gather_recvcnt
      glb_field( comm%gather_recvmap(jnode) ) = recvbuffer(jnode)
    end do

  end subroutine gather_real8_rank1


  subroutine gather_real8_rank2(comm,loc_field,glb_field)
    class(Comm_type), intent(inout) :: comm
    real*8, dimension(:,:), intent(in) :: loc_field
    real*8, dimension(:,:), allocatable, intent(inout) :: glb_field
    real*8 :: sendbuffer(comm%gather_sendcnt)
    real*8 :: recvbuffer(comm%gather_recvcnt)
    integer :: jnode, ncols, jcol

    ncols = size(loc_field,2)

    if( .not. allocated( glb_field) ) then
      allocate( glb_field( comm%glb_field_size(), ncols ) )
    end if

    do jcol=1,ncols

      ! Pack
      do jnode=1,comm%gather_sendcnt
        sendbuffer(jnode) = loc_field(comm%gather_sendmap(jnode),jcol)
      end do

      ! Communicate
      call MPI_GATHERV( sendbuffer,comm%gather_sendcnt,MPI_DOUBLE_PRECISION, &
                      & recvbuffer,comm%gather_recvcounts,comm%gather_recvdispls,MPI_DOUBLE_PRECISION, &
                      & comm%gather_root,MPI_COMM_WORLD,ierr )
    
      ! Unpack
      do jnode=1,comm%gather_recvcnt
         glb_field( comm%gather_recvmap(jnode), jcol ) = recvbuffer(jnode)
      end do
    end do

  end subroutine gather_real8_rank2

  function glb_field_size(comm) result(rows)
    class(Comm_type), intent(inout) :: comm
    integer :: rows
    rows = comm%gather_recvcnt
  end function glb_field_size

end module parallel_module
