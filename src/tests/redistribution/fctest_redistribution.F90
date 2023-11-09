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

function do_redistribute(fspace_1, fspace_2) result(field_2)
use atlas_module
use atlas_redistribution_module

implicit none

type(atlas_FunctionSpace), intent(in) :: fspace_1, fspace_2
type(atlas_Field) :: field_2

type(atlas_FunctionSpace) :: fspace_hlp
type(atlas_Redistribution) :: redist, redist_hlp
type(atlas_Field) :: field_1 !, field_2
real(c_float), pointer :: field_1v(:), field_2v(:)

redist = atlas_Redistribution(fspace_1, fspace_2)

field_1 = fspace_1%create_field(atlas_real(c_float))
field_2 = fspace_2%create_field(atlas_real(c_float))
call field_1%data(field_1v)
field_1v = 1._c_float
call field_2%data(field_2v)
field_2v = 2._c_float

call redist%execute(field_1, field_2)

call field_2%data(field_2v)
call field_2%halo_exchange()

! check access to source and target function spaces
redist_hlp = atlas_Redistribution(redist%c_ptr())
fspace_hlp = redist%source()
fspace_hlp = redist%target()

call field_1%final()
call redist_hlp%final()
call redist%final()
end function do_redistribute

end module

! -----------------------------------------------------------------------------

TESTSUITE(fcta_Redistribution)
!TESTSUITE_WITH_FIXTURE(fcta_Redistribution,fcta_Redistribution_fxt)

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

TEST( test_redistribution )
use atlas_module
use fcta_Redistribution_fxt
use atlas_redistribution_module

implicit none

type(atlas_Grid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_FunctionSpace)  :: fspace_1, fspace_2
type(atlas_Field) :: field
real(c_float) , pointer :: field_v(:)

grid = atlas_Grid("O8")

fspace_1 = atlas_functionspace_StructuredColumns(grid, atlas_Partitioner("equal_regions"))
fspace_2 = atlas_functionspace_StructuredColumns(grid, atlas_Partitioner("regular_bands"))
field = do_redistribute(fspace_1, fspace_2)
call field%data(field_v)
FCTEST_CHECK(maxval(field_v) == 1.)

meshgenerator = atlas_MeshGenerator()
mesh = meshgenerator%generate(grid, atlas_Partitioner("equal_regions"))
fspace_1 = atlas_functionspace_NodeColumns(mesh)
mesh = meshgenerator%generate(grid, atlas_Partitioner("regular_bands"))
fspace_2 = atlas_functionspace_NodeColumns(mesh)
field = do_redistribute(fspace_1, fspace_2)
call field%data(field_v)
FCTEST_CHECK(maxval(field_v) == 1.)

fspace_2 = atlas_functionspace_EdgeColumns(mesh)
mesh = meshgenerator%generate(grid, atlas_Partitioner("equal_regions"))
fspace_1 = atlas_functionspace_EdgeColumns(mesh)
field = do_redistribute(fspace_1, fspace_2)
call field%data(field_v)
FCTEST_CHECK(maxval(field_v) == 1.)

call fspace_2%final()
call fspace_1%final()
call grid%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
