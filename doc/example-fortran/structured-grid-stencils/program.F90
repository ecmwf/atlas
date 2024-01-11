! (C) Copyright 2023 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------------------------
! Example program on how to access and use StructuredColumns functionspace
! Also includes computation of stencils around a given coordinate
! ------------------------------------------------------------------------------------------------

subroutine run

use, intrinsic :: iso_fortran_env, only : real64, real32
use atlas_module, only :   &
  & atlas_StructuredGrid  ,&
  & atlas_functionspace_StructuredColumns ,&
  & atlas_Field, &
  & ATLAS_KIND_IDX, &
  & atlas_StructuredGrid_ComputeNorth ,&
  & atlas_StructuredGrid_ComputeWest ,&
  & atlas_StructuredGrid_ComputeStencil ,&
  & atlas_StructuredGrid_Stencil, &
  & atlas_real

use fckit_mpi_module, only : fckit_mpi

! ------------------------------------------------
implicit none

character(len=:), allocatable :: gridname
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_StructuredColumns) :: functionspace
real(real64) :: zlon, zlat
integer(ATLAS_KIND_IDX) :: jlon,jlat, jgp
real(real64), pointer :: xy(:,:)
integer :: log_rank

! ------------------------------------------------

log_rank = 1 ! First partition has rank=0 !!!
if (fckit_mpi%size() == 1) then
  log_rank = 0
endif

gridname = "O320"

! Create grid and write some properties
grid          = atlas_StructuredGrid(gridname)

if (fckit_mpi%rank() == log_rank) then
    write(0,*) "grid%name()  : ", grid%name()
    write(0,*) "grid%size()  : ", grid%size()
    write(0,*) "grid%ny()    : ", grid%ny()
    write(0,*) "grid%nxmax() : ", grid%nxmax()
endif

! Create functionspace of grid, distributed across MPI tasks with default partitioner, and with halo of 2 wide
!   The default partitioner is "equal_regions"
functionspace = atlas_functionspace_StructuredColumns(grid, halo=2)

! Latitude bounds without and with halos of this partition
if (fckit_mpi%rank() == log_rank) then
    write(0,*) "functionspace%size_owned() : ", functionspace%size_owned()
    write(0,*) "functionspace%size()       : ", functionspace%size()
    write(0,'(A,I0,A,I0)') " halo points between index ", functionspace%size_owned()+1 , " and " , functionspace%size()
    write(0,*) "functionspace%j_begin(),      functionspace%j_end()      : ", &
      & functionspace%j_begin(), functionspace%j_end()
    write(0,*) "functionspace%j_begin_halo(), functionspace%j_end_halo() : ", &
      & functionspace%j_begin_halo(), functionspace%j_end_halo()
endif

! Longitude bounds without and with halos of this partition for the first latitude
jlat = functionspace%j_begin()
if (fckit_mpi%rank() == log_rank) then
  write(0,*) "functionspace%i_begin(jlat),      functionspace%i_end(jlat)      : ", functionspace%i_begin(jlat),     &
   & functionspace%i_end(jlat)
  write(0,*) "functionspace%i_begin_halo(jlat), functionspace%i_end_halo(jlat) : ", functionspace%i_begin_halo(jlat),&
   & functionspace%i_end_halo(jlat)
endif

! Access i,j indices of grid points
block
  type(atlas_Field) :: field_index_i, field_index_j
  integer(ATLAS_KIND_IDX), pointer :: index_i(:), index_j(:)
  field_index_i = functionspace%index_i()
  field_index_j = functionspace%index_j()
  call field_index_i%data(index_i)
  call field_index_j%data(index_j)
  if (fckit_mpi%rank() == log_rank) then
    write(0,*) "i,j of first partition point :", index_i(1), index_j(1)
  endif
  call field_index_i%final()
  call field_index_j%final()
end block

! Access to xy coordinates
block
  type(atlas_Field) :: field_xy
  field_xy = functionspace%xy()
  call field_xy%data(xy)
  call field_xy%final()
end block

! Creating a horizontal field and perform halo exchange
block
  type(atlas_Field) :: field
  real(real32), pointer :: view(:)
  field = functionspace%create_field(name="myfield", kind=atlas_real(real32))
  call field%data(view)
  if (fckit_mpi%rank() == log_rank) then
    write(0,*) "shape( horizontal field )", shape(view)
  endif
  view(1:functionspace%size_owned()) = 1.
  view(functionspace%size_owned()+1:functionspace%size()) = 0
  call field%set_dirty() ! This marks that halos are in "dirty" state
  call field%halo_exchange() ! Only exchanges halos when halos are marked as "dirty"
  if (fckit_mpi%rank() == log_rank) then
    write(0,*) "halo exhange success : ", all(view == 1.)
  endif
  call field%final()
end block

! Creating a horizontal/vertical field
block
  type(atlas_Field) :: field
  real(real32), pointer :: view(:,:)
  field = functionspace%create_field(name="myfield", kind=atlas_real(real32), levels=10)
  call field%data(view)
  if (fckit_mpi%rank() == log_rank) then
    write(0,*) "shape( horizontal/vertical field )", shape(view)
  endif
  call field%final()
end block

! Set a coordinate somewhere south-east of the first grid point, but still within this partition
jgp = 1
zlon = xy(1,jgp) + 0.1
zlat = xy(2,jgp) - 0.1

! Compute nearest points to the north-west of a coordinate
block
  type(atlas_StructuredGrid_ComputeNorth) :: compute_north
  type(atlas_StructuredGrid_ComputeWest) :: compute_west
  compute_north = atlas_StructuredGrid_ComputeNorth(grid, halo=2)
  compute_west  = atlas_StructuredGrid_ComputeWest(grid,  halo=2)
  
  jlat = compute_north%execute(zlat)
  jlon = compute_west%execute(zlon, jlat)
  jgp = functionspace%index(jlon,jlat) ! gridpoint index in 1-dimensional array
  
  if (fckit_mpi%rank() == log_rank) then
  write(0,'(A,F6.2,A,I0)') "compute_north%execute(y=",zlat,") : ", jlat
  write(0,'(A,F6.2,A,I0,A,I0)') "compute_west%execute(x=",zlon,", j=", jlat,") : ", jlon
  write(0,'(A,I0,A,F6.2,A,F6.2,A)') "gridpoint north-west: jpg=",jgp,"    (lon,lat)=(", xy(1,jgp), ",", xy(2,jgp), ")"
  endif
  call compute_west%final()
  call compute_north%final()
end block

! Stencil computer of 4x4 stencil around a coordinate
block
  type(atlas_StructuredGrid_ComputeStencil) :: compute_stencil
  type(atlas_StructuredGrid_Stencil) :: stencil
  call compute_stencil%setup(grid, stencil_width=4)
  call compute_stencil%execute(zlon,zlat,stencil)
  if (fckit_mpi%rank() == log_rank) then
    write(0,'(A,A,dt)') 'stencil:', new_line('a'), stencil
  endif
  if (fckit_mpi%rank() == log_rank) then
    write(0,'(A,A,dt)') 'stencil gridpoints:'
    do jlat=1,stencil%width
      do jlon=1,stencil%width
        write(0,'(I8)',advance='no') functionspace%index(stencil%i(jlon,jlat),stencil%j(jlat))
      enddo
      write(0,'(A)') ""
    enddo

    write(0,'(A,A,dt)') 'stencil coordinates:'
    do jlat=1,stencil%width
      jgp = functionspace%index(stencil%i(1,jlat),stencil%j(jlat))
      write(0,'(A, F6.2,A)', advance='no') " lat : ", xy(2,jgp), "  ---  lon : "
      do jlon=1,stencil%width
        jgp = functionspace%index(stencil%i(jlon,jlat),stencil%j(jlat))
        write(0,'(F8.2)',advance='no') xy(1,jgp)
      enddo
      write(0,'(A)') ""
    enddo

  endif
  call compute_stencil%final()
end block

call functionspace%final()
call grid%final()
end subroutine

! ------------------------------------------------------------------------------------------------

program main
use atlas_module, only : atlas_library
implicit none
call atlas_library%initialize()
call run()
call atlas_library%finalize()
end program