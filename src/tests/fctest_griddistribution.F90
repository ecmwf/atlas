! (C) Copyright 1996-2016 ECMWF.
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

! -----------------------------------------------------------------------------

TEST( test_griddist )
  use atlas_module
  implicit none
  type(atlas_grid_Structured) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_GridDistribution) :: griddistribution

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
  meshgenerator = atlas_meshgenerator_Structured()
  mesh = meshgenerator%generate(grid,griddistribution)
  call griddistribution%final()

  call atlas_write_gmsh(mesh,"testf3.msh")

  deallocate(part)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

