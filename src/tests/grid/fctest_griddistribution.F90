! (C) Copyright 1996-2017 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

TESTSUITE(fctest_GridDist)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  use atlas_module
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use atlas_module
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! ----------------------------------------------------------------------------

! T E S T( test_reducedgaussian )
!   use atlas_module
!   use, intrinsic :: iso_c_binding
!   implicit none
!   type(atlas_grid_Structured) :: N640
!   type(atlas_grid_ReducedGaussian) :: custom
!   integer(c_long), pointer :: pl(:)
!
!   N640 = atlas_grid_Structured("N640")
!   FCTEST_CHECK_EQUAL(N640%npts(),2140702_c_long)
!   pl => N640%pl()
!
!   custom = atlas_grid_ReducedGaussian( N640%N(), pl )
!   FCTEST_CHECK_EQUAL(N640%npts(),custom%npts() )
!
!   call N640%final()
!   call custom%final()
!
! E N D _ T E S T

! -----------------------------------------------------------------------------

TEST( test_griddist )
  use atlas_module
  implicit none
  type(atlas_grid_Structured) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_Output) :: gmsh
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_GridDistribution) :: griddistribution
  character(len=1024) :: msg

  integer, allocatable :: part(:)
  integer :: jnode

  grid = atlas_grid_Structured("O16")
  !grid = atlas_grid_Structured("ll.128x64")
  !grid = atlas_grid_ShiftedLonLat(128,64)

  allocate( part(grid%npts()) )
  do jnode=1,grid%npts()/3
    part(jnode) = 1
  enddo
  do jnode=grid%npts()/3+1,grid%npts()
    part(jnode) = 1
  enddo

  griddistribution = atlas_GridDistribution(part, part0=1)

  write(msg,*) "owners fort",griddistribution%owners()
  call atlas_log%info(msg)

  meshgenerator = atlas_meshgenerator_Structured()
  mesh = meshgenerator%generate(grid,griddistribution)
  call griddistribution%final()

  gmsh = atlas_output_Gmsh("testf3.msh")
  call gmsh%write(mesh)

  deallocate(part)

  call mesh%final()
  call gmsh%final()
  call grid%final()
  call meshgenerator%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

