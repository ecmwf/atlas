! (C) Copyright 2013 ECMWF.
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
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use atlas_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! ----------------------------------------------------------------------------

! T E S T( test_reducedgaussian )
!   use atlas_module
!   use, intrinsic :: iso_c_binding
!   implicit none
!   type(atlas_StructuredGrid) :: N640
!   type(atlas_ReducedGaussianGrid) :: custom
!   integer(c_long), pointer :: pl(:)
!
!   N640 = atlas_StructuredGrid("N640")
!   FCTEST_CHECK_EQUAL(N640%size(),2140702_c_long)
!   pl => N640%pl()
!
!   custom = atlas_ReducedGaussianGrid( N640%N(), pl )
!   FCTEST_CHECK_EQUAL(N640%size(),custom%size() )
!
!   call N640%final()
!   call custom%final()
!
! E N D _ T E S T

! -----------------------------------------------------------------------------

TEST( test_griddist )
  use atlas_module
  implicit none
  type(atlas_StructuredGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_Output) :: gmsh
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_GridDistribution) :: griddistribution

  integer, allocatable :: part(:)
  integer(ATLAS_KIND_IDX) :: jnode

  grid = atlas_StructuredGrid("O16")

  allocate( part(grid%size()) )
  do jnode=1,grid%size()/3
    part(jnode) = 1
  enddo
  do jnode=grid%size()/3+1,grid%size()
    part(jnode) = 1
  enddo

  griddistribution = atlas_GridDistribution(part, part0=1)

  FCTEST_CHECK_EQUAL( griddistribution%nb_partitions(), 1 )
  FCTEST_CHECK_EQUAL( griddistribution%nb_pts(), [ int(grid%size(),ATLAS_KIND_IDX) ] )


  FCTEST_CHECK_EQUAL( grid%owners(), 1 )

  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid,griddistribution)
  FCTEST_CHECK_EQUAL( mesh%owners(), 1 )

  call griddistribution%final()

  FCTEST_CHECK_EQUAL( grid%owners(), 2 )

  gmsh = atlas_output_Gmsh("testf3.msh")
  call gmsh%write(mesh)

  deallocate(part)

  call mesh%final()
  FCTEST_CHECK_EQUAL( grid%owners(), 1 )
  call gmsh%final()
  call grid%final()
  call meshgenerator%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

