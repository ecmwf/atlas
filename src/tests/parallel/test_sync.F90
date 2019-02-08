! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ===================================================================
! test_sync program
! ---------------------
! This program tests MPI_ALLTOALLV communication,
! used for synchronisation
! ===================================================================
program test_sync

  use common_module
  use parallel_module
  use datastruct

  implicit none

  integer, allocatable  :: nb_nodes(:)
  integer, allocatable  :: proc(:)
  integer, allocatable  :: glb_idx(:)
  integer, allocatable  :: master_glb_idx(:)
  real(kind=jprs), allocatable, target  :: field(:), vectorfield(:,:)
  real(kind=jprs), pointer :: vectorfield_ptr(:,:)
  real(kind=jprs), allocatable  :: glb_field(:)
  integer :: bounds(2)
  integer :: length
  integer :: parallel_bound = 1
  type(Comm_type) :: comm

  type(HaloExchange_type) :: halo_exchange
  type(FunctionSpace_type) :: function_space

  call set_log_level(LOG_LEVEL_INFO)

  call parallel_init() ! MPI_INIT etc

  halo_exchange = new_HaloExchange()

  ! Create distribution with overlap regions

  if (nproc .eq. 3) then
    allocate( nb_nodes(3) )

    nb_nodes = [ 5, 6, 7 ]
    length = nb_nodes(myproc+1)
    bounds = [ -1, length ]
    allocate( proc(length) )
    allocate( glb_idx(length) )
    allocate( master_glb_idx(length) )
    allocate( field(length) )
    allocate( vectorfield(2,length) )

    select case (myproc)
      case (0)
        proc = [2,0,0,0,1]
        glb_idx = [9,1,2,3,4]
        field = [-1,1,2,3,-1]
      case (1)
        proc = [0,1,1,1,2,2]
        glb_idx = [3,4,5,6,7,8]
        field = [-1,4,5,6,-1,-1]
      case (2)
        proc = [1,1,2,2,2,0,0]
        glb_idx = [5,6,7,8,9,1,2]
        field = [-1,-1,7,8,9,-1,-1]
    end select
    master_glb_idx(:) = glb_idx(:)
    vectorfield(1,:) = field*10
    vectorfield(2,:) = field*100

    vectorfield_ptr => vectorfield
    ! Setup a communicator for synchronisation
    call comm%setup(proc,glb_idx)

    ! We can setup function_space for halo_exchange
    function_space = new_FunctionSpace("nodes","shape_func",length)
    call function_space%parallelise( proc, glb_idx, master_glb_idx )

    ! Or we can setup a custom halo_exchange object
    !call halo_exchange%setup(proc, glb_idx, bounds, parallel_bound)

    ! Update the field values whose proc is not myproc
    ! This is the older fortran implementation
    !call comm%synchronise(field)
    !call comm%synchronise(vectorfield)

    ! Halo exchange through the function_space (We don't need to know nbvars)
    call function_space%halo_exchange(field)
    call function_space%halo_exchange(vectorfield)
    !call function_space%halo_exchange(vectorfield(:))

    !call function_space%halo_exchange(vectorfield(:,2))

    ! Halo exchange through the custom halo_exchange object
    !call halo_exchange%execute(field,1)
    !call halo_exchange%execute(vectorfield,2)
    !call halo_exchange%execute(vectorfield(:,1),1)
    !call halo_exchange%execute(vectorfield(:,2),1)

    ! Verify that update happened correctly
    write(log_str,*) myproc, ": field            = ", field; call log_info()
    write(log_str,*) myproc, ": vectorfield(1,:) = ", vectorfield(1,:); call log_info()
    write(log_str,*) myproc, ": vectorfield(2,:) = ", vectorfield(2,:); call log_info()

    !call comm%gather(field, glb_field)

    !call set_log_proc(0)
    !write(log_str,*) myproc, ": glb_field        = ", glb_field; call log_info()

  else
    call set_log_proc(0)
    call log_error("ERROR: Run this program with 3 tasks")
  end if

  call delete(halo_exchange)

  call parallel_finalise() ! MPI_FINALIZE

end program test_sync

! OUTPUT:
!
! 0: field =  90    10    20    30   40
! 1: field =  30    40    50    60   70
! 2: field =  60    70    80    90   10
