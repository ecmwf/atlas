! ===================================================================
! test_sync program
! ---------------------
! This program tests MPI_ALLTOALLV communication, 
! used for synchronisation
! ===================================================================
program test_sync

  use common_module
  use parallel_module

  implicit none

  integer, allocatable  :: proc(:)
  integer, allocatable  :: glb_idx(:)
  real(kind=jprb), allocatable  :: field(:), vectorfield(:,:)
  real(kind=jprb), allocatable  :: glb_field(:)
  integer :: length 
  type(Comm_type) :: comm

  call set_log_level(LOG_LEVEL_INFO)
  call set_log_proc(0)

  call parallel_init() ! MPI_INIT etc

  ! Create distribution with overlap regions
  length = 5
  allocate( proc(length) )
  allocate( glb_idx(length) )
  allocate( field(length) )
  allocate( vectorfield(length,2) )

  if (nproc .eq. 3) then
  
    select case (myproc)
      case (0)
        proc = [2,0,0,0,1]
        glb_idx = [9,1,2,3,4]
        field = [-1,10,20,30,-1]
      case (1)
        proc = [0,1,1,1,2]
        glb_idx = [3,4,5,6,7]
        field = [-1,40,50,60,-1]
      case (2)
        proc = [1,2,2,2,0]
        glb_idx = [6,7,8,9,1]
        field = [-1,70,80,90,-1]
    end select
    vectorfield(:,1) = field
    vectorfield(:,2) = field

    ! Setup a communicator for synchronisation
    call comm%setup(proc,glb_idx)

    ! Update the field values whose proc is not myproc
    call comm%synchronise(field)
    call comm%synchronise(vectorfield)
  
    ! Verify that update happened correctly
    write(log_str,*) myproc, ": field = ", field; call log_info()

    call comm%gather(field, glb_field)

    write(log_str,*) myproc, ": glb_field = ", glb_field; call log_info()

  else
    call set_log_proc(0)
    call log_error("ERROR: Run this program with 3 tasks")
  end if

  call parallel_finalise() ! MPI_FINALIZE

end program test_sync

! OUTPUT:
!
! 0: field =  90    10    20    30   40
! 1: field =  30    40    50    60   70
! 2: field =  60    70    80    90   10
