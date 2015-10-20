! (C) Copyright 1996-2015 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_nabla_EdgeBasedFiniteVolume_Fixture
use atlas_module
use iso_c_binding
implicit none

contains

end module fctest_atlas_nabla_EdgeBasedFiniteVolume_Fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_nabla_EdgeBasedFiniteVolume,fctest_atlas_nabla_EdgeBasedFiniteVolume_Fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_fvm )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_functionspace_EdgeBasedFiniteVolume) :: fvm

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fvm  = atlas_functionspace_EdgeBasedFiniteVolume(mesh)

call fvm%final()
call mesh%final()
call grid%final()

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

