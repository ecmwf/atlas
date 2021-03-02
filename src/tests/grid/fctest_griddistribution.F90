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

! -----------------------------------------------------------------------------

TEST( test_griddist )
#if 1
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

  deallocate(part)

  FCTEST_CHECK_EQUAL( grid%owners(), 2 )

  gmsh = atlas_output_Gmsh("testf3.msh")
  call gmsh%write(mesh)

  do jnode=1,grid%size()
    FCTEST_CHECK_EQUAL( griddistribution%partition(jnode), 0 )
  enddo

  call griddistribution%final()

  call mesh%final()
  FCTEST_CHECK_EQUAL( grid%owners(), 1 )
  call gmsh%final()
  call grid%final()
  call meshgenerator%final()
#endif
END_TEST

! -----------------------------------------------------------------------------

TEST( test_griddist_from_partitioner )
  use atlas_module
  implicit none
  type(atlas_StructuredGrid) :: grid
  type(atlas_Partitioner) :: partitioner
  type(atlas_GridDistribution) :: griddistribution

  integer(ATLAS_KIND_IDX) :: jnode

  grid = atlas_StructuredGrid("O16")

  griddistribution = atlas_GridDistribution(grid,atlas_Partitioner("serial"))

  FCTEST_CHECK_EQUAL( griddistribution%nb_partitions(), 1 )
  FCTEST_CHECK_EQUAL( griddistribution%nb_pts(), [ int(grid%size(),ATLAS_KIND_IDX) ] )

  do jnode=1,grid%size()
    FCTEST_CHECK_EQUAL( griddistribution%partition(jnode), 0 )
  enddo

  call griddistribution%final()
  call grid%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

