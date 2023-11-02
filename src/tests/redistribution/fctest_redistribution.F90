! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
!
! @author Willem Deconinck
! @author Slavko Brdar

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_Redistribution_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg
contains

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_Redistribution,fcta_Redistribution_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

TEST( test_resitribution )
use atlas_module
use atlas_redistribution_module

implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_StructuredColumns) :: fspace1, fspace2
type(atlas_Redistribution) :: redist
type(c_ptr) :: cptr
grid = atlas_StructuredGrid("O16")
fspace1 = atlas_functionspace_StructuredColumns(grid, atlas_Partitioner("equal_regions"), halo=2)
fspace2 = atlas_functionspace_StructuredColumns(grid, atlas_Partitioner("regular_bands"))

redist = atlas_Redistribution(fspace1, fspace2)

call redist%final()
call fspace2%final()
call fspace1%final()
call grid%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

