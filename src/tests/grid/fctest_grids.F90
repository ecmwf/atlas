! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

TESTSUITE(fctest_Grid)

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

TEST( test_ij2gidx )
  use atlas_module
  implicit none
  type(atlas_StructuredGrid) :: grid
  integer(ATLAS_KIND_IDX) :: i, j, i1, j1
  integer(ATLAS_KIND_IDX) :: jglo, jglo1
  grid = atlas_StructuredGrid ("N16")

  jglo = 1
  do j = 1, grid%ny ()
  do i = 1, grid%nx (j)
    call grid%index2ij (jglo, i1, j1)
    jglo1 = grid%index (i, j)
    FCTEST_CHECK_EQUAL( jglo, jglo1 )
    FCTEST_CHECK_EQUAL( i, i1 )
    FCTEST_CHECK_EQUAL( j, j1 )
    jglo = jglo + 1
  enddo
  enddo

  call grid%final()
  
END_TEST

TEST( test_spec )
  use atlas_module
  implicit none
  type(atlas_Grid) :: grid
  type(atlas_Config) :: spec
  character(:), allocatable :: json_sorted, json_ordered, json

  grid = atlas_Grid ("O32")
  spec = grid%spec()

  FCTEST_CHECK_EQUAL( spec%owners(), 1 )
  json_sorted  = '{"domain":{"type":"global"},"name":"O32","projection":{"type":"lonlat"}}'
  json_ordered = '{"name":"O32","domain":{"type":"global"},"projection":{"type":"lonlat"}}'
  json = spec%json()
  FCTEST_CHECK( json == json_sorted .or. json == json_ordered )

  grid = atlas_RegularGaussianGrid(8)
  grid = atlas_RegularLonLatGrid(8,8)
  grid = atlas_ShiftedLonLatGrid(8,8)
  grid = atlas_ShiftedLonGrid(8,8)
  grid = atlas_ShiftedLatGrid(8,8)

  call spec%final()
  call grid%final()

END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

