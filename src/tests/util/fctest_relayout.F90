! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "fckit/fctest.h"

!> @file fctest_relayout.F90
!> @brief Unit tests for the Fortran atlas_Relayout_module interfaces.
!>
!> The tests fill structured and blocked fields with checksum-compatible values, call the
!> Fortran relayout generics, and verify the result with StructuredColumns and
!> BlockStructuredColumns checksums.

! -----------------------------------------------------------------------------

module fctest_Relayout_fxt
use atlas_module, only : atlas_library, atlas_StructuredGrid, atlas_Field, atlas_FieldSet, atlas_real, atlas_kind_idx
use atlas_functionspace_BlockStructuredColumns_module
use atlas_functionspace_StructuredColumns_module
use atlas_Relayout_module, only : copy_blocked_to_blocked, copy_blocked_to_nonblocked, copy_nonblocked_to_blocked
use, intrinsic :: iso_c_binding
implicit none

contains

!> @brief Create a double precision blocked field with optional levels and variables.
!>
!> @param functionspace BlockStructuredColumns functionspace used to create the field.
!> @param name Field name.
!> @param levels Number of levels.  Must be non-negative.
!> @param variables Number of variables.  Must be non-negative.  If variables is positive,
!>        levels must also describe the middle vertical dimension used by the relayout tests.
function create_blocked_field(functionspace, name, levels, variables) result(field)
  type(atlas_functionspace_BlockStructuredColumns), intent(in) :: functionspace
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  integer, intent(in) :: variables
  type(atlas_Field) :: field

  if (variables > 0) then
    field = functionspace%create_field(name=name, kind=atlas_real(c_double), levels=levels, variables=variables)
  else if (levels > 0) then
    field = functionspace%create_field(name=name, kind=atlas_real(c_double), levels=levels)
  else
    field = functionspace%create_field(name=name, kind=atlas_real(c_double))
  endif
end function create_blocked_field

!> @brief Create a double precision structured field with optional levels and variables.
!>
!> @param functionspace StructuredColumns functionspace used to create the field.
!> @param name Field name.
!> @param levels Number of levels.  Must be non-negative.
!> @param variables Number of variables.  Must be non-negative.  If variables is positive,
!>        levels must also describe the middle vertical dimension used by the relayout tests.
function create_structured_field(functionspace, name, levels, variables) result(field)
  type(atlas_functionspace_StructuredColumns), intent(in) :: functionspace
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  integer, intent(in) :: variables
  type(atlas_Field) :: field

  if (variables > 0) then
    field = functionspace%create_field(name=name, kind=atlas_real(c_double), levels=levels, variables=variables)
  else if (levels > 0) then
    field = functionspace%create_field(name=name, kind=atlas_real(c_double), levels=levels)
  else
    field = functionspace%create_field(name=name, kind=atlas_real(c_double))
  endif
end function create_structured_field

!> @brief Fill a structured field in checksum order.
!>
!> @param field Structured field with rank 1, 2, or 3 and real(c_double) data.
!> @param first_value First value written to the field.
subroutine fill_structured_field(field, first_value)
  type(atlas_Field), intent(inout) :: field
  real(c_double), intent(in) :: first_value
  real(c_double), pointer :: values_1d(:)
  real(c_double), pointer :: values_2d(:,:)
  real(c_double), pointer :: values_3d(:,:,:)
  real(c_double) :: next_value
  integer :: jpoint, jlev, jvar, jentry

  next_value = first_value
  if (field%rank() == 3) then
    call field%data(values_3d)
    do jpoint = 1, size(values_3d, 1)
      do jlev = 1, size(values_3d, 2)
        do jvar = 1, size(values_3d, 3)
          values_3d(jpoint, jlev, jvar) = next_value
          next_value = next_value + 1._c_double
        enddo
      enddo
    enddo
  else if (field%rank() == 2) then
    call field%data(values_2d)
    do jpoint = 1, size(values_2d, 1)
      do jentry = 1, size(values_2d, 2)
        values_2d(jpoint, jentry) = next_value
        next_value = next_value + 1._c_double
      enddo
    enddo
  else
    call field%data(values_1d)
    do jpoint = 1, size(values_1d, 1)
      values_1d(jpoint) = next_value
      next_value = next_value + 1._c_double
    enddo
  endif
end subroutine fill_structured_field

!> @brief Fill a blocked field in an order that gives the same checksum as a structured field.
!>
!> @param functionspace BlockStructuredColumns functionspace that owns `field`.
!> @param field Blocked field with rank 2, 3, or 4 and real(c_double) data.
!> @param first_value First value written to the field.
!>
!> @pre `field` must have been created by `functionspace` so that `block_size(jblk)` matches
!>      the field block dimension.
subroutine fill_blocked_field(functionspace, field, first_value)
  type(atlas_functionspace_BlockStructuredColumns), intent(in) :: functionspace
  type(atlas_Field), intent(inout) :: field
  real(c_double), intent(in) :: first_value
  real(c_double), pointer :: values_2d(:,:)
  real(c_double), pointer :: values_3d(:,:,:)
  real(c_double), pointer :: values_4d(:,:,:,:)
  real(c_double) :: next_value
  integer(ATLAS_KIND_IDX) :: jblk, jlane, jlev, jvar, jentry

  next_value = first_value
  if (field%rank() == 4) then
    call field%data(values_4d)
    do jblk = 1, size(values_4d, 4)
      do jlane = 1, functionspace%block_size(jblk)
        do jlev = 1, size(values_4d, 3)
          do jvar = 1, size(values_4d, 2)
            values_4d(jlane, jvar, jlev, jblk) = next_value
            next_value = next_value + 1._c_double
          enddo
        enddo
      enddo
    enddo
  else if (field%rank() == 3) then
    call field%data(values_3d)
    do jblk = 1, size(values_3d, 3)
      do jlane = 1, functionspace%block_size(jblk)
        do jentry = 1, size(values_3d, 2)
          values_3d(jlane, jentry, jblk) = next_value
          next_value = next_value + 1._c_double
        enddo
      enddo
    enddo
  else
    call field%data(values_2d)
    do jblk = 1, size(values_2d, 2)
      do jlane = 1, functionspace%block_size(jblk)
        values_2d(jlane, jblk) = next_value
        next_value = next_value + 1._c_double
      enddo
    enddo
  endif
end subroutine fill_blocked_field

!> @brief Check the Fortran blocked-to-blocked relayout wrapper for one field shape.
!>
!> @param source_fs Source blocked functionspace.
!> @param target_fs Target blocked functionspace.  It may use a different nproma.
!> @param levels Number of levels used to create the test field.
!> @param variables Number of variables used to create the test field.
!> @param first_value First source value.
function check_blocked_to_blocked(source_fs, target_fs, levels, variables, first_value)
  logical :: check_blocked_to_blocked
  type(atlas_functionspace_BlockStructuredColumns), intent(in) :: source_fs
  type(atlas_functionspace_BlockStructuredColumns), intent(in) :: target_fs
  integer, intent(in) :: levels
  integer, intent(in) :: variables
  real(c_double), intent(in) :: first_value
  type(atlas_Field) :: source, target

  source = create_blocked_field(source_fs, "source", levels, variables)
  target = create_blocked_field(target_fs, "target", levels, variables)

  call fill_blocked_field(source_fs, source, first_value)
  call copy_blocked_to_blocked(source, target)

  check_blocked_to_blocked = (source_fs%checksum(source) == target_fs%checksum(target))
  !FCTEST_CHECK_EQUAL(source_fs%checksum(source), target_fs%checksum(target))
end function check_blocked_to_blocked

!> @brief Check the Fortran blocked-to-nonblocked relayout wrapper for one field shape.
!>
!> @param structured_fs Target structured functionspace.
!> @param block_fs Source blocked functionspace.
!> @param levels Number of levels used to create the test field.
!> @param variables Number of variables used to create the test field.
!> @param first_value First source value.
function check_blocked_to_nonblocked(structured_fs, block_fs, levels, variables, first_value)
  logical :: check_blocked_to_nonblocked
  type(atlas_functionspace_StructuredColumns), intent(in) :: structured_fs
  type(atlas_functionspace_BlockStructuredColumns), intent(in) :: block_fs
  integer, intent(in) :: levels
  integer, intent(in) :: variables
  real(c_double), intent(in) :: first_value
  type(atlas_Field) :: source, target

  source = create_blocked_field(block_fs, "source", levels, variables)
  target = create_structured_field(structured_fs, "target", levels, variables)

  call fill_blocked_field(block_fs, source, first_value)
  call copy_blocked_to_nonblocked(source, target)

  check_blocked_to_nonblocked = (block_fs%checksum(source) == structured_fs%checksum(target))
  !FCTEST_CHECK_EQUAL(block_fs%checksum(source), structured_fs%checksum(target))
end function check_blocked_to_nonblocked

!> @brief Check the Fortran nonblocked-to-blocked relayout wrapper for one field shape.
!>
!> @param structured_fs Source structured functionspace.
!> @param block_fs Target blocked functionspace.
!> @param levels Number of levels used to create the test field.
!> @param variables Number of variables used to create the test field.
!> @param first_value First source value.
function check_nonblocked_to_blocked(structured_fs, block_fs, levels, variables, first_value)
  logical :: check_nonblocked_to_blocked
  type(atlas_functionspace_StructuredColumns), intent(in) :: structured_fs
  type(atlas_functionspace_BlockStructuredColumns), intent(in) :: block_fs
  integer, intent(in) :: levels
  integer, intent(in) :: variables
  real(c_double), intent(in) :: first_value
  type(atlas_Field) :: source, target

  source = create_structured_field(structured_fs, "source", levels, variables)
  target = create_blocked_field(block_fs, "target", levels, variables)

  call fill_structured_field(source, first_value)
  call copy_nonblocked_to_blocked(source, target)

  check_nonblocked_to_blocked = (structured_fs%checksum(source) == block_fs%checksum(target))
  !FCTEST_CHECK_EQUAL(structured_fs%checksum(source), block_fs%checksum(target))
end function check_nonblocked_to_blocked

end module fctest_Relayout_fxt

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_Relayout, fctest_Relayout_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST(test_relayout_fields)
implicit none
  type(atlas_StructuredGrid) :: grid
  type(atlas_functionspace_StructuredColumns) :: structured_fs
  type(atlas_functionspace_BlockStructuredColumns) :: source_block_fs
  type(atlas_functionspace_BlockStructuredColumns) :: target_block_fs

  grid = atlas_StructuredGrid("O8")
  structured_fs = atlas_functionspace_StructuredColumns(grid, halo=0)
  source_block_fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=5)
  target_block_fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=7)

  FCTEST_CHECK(check_blocked_to_blocked(source_block_fs, target_block_fs, 0, 0, 1._c_double))
  FCTEST_CHECK(check_blocked_to_blocked(source_block_fs, target_block_fs, 4, 0, 1000._c_double))
  FCTEST_CHECK(check_blocked_to_blocked(source_block_fs, target_block_fs, 4, 3, 2000._c_double))

  FCTEST_CHECK(check_blocked_to_nonblocked(structured_fs, source_block_fs, 0, 0, 3000._c_double))
  FCTEST_CHECK(check_blocked_to_nonblocked(structured_fs, source_block_fs, 4, 0, 4000._c_double))
  FCTEST_CHECK(check_blocked_to_nonblocked(structured_fs, source_block_fs, 4, 3, 5000._c_double))

  FCTEST_CHECK(check_nonblocked_to_blocked(structured_fs, source_block_fs, 0, 0, 6000._c_double))
  FCTEST_CHECK(check_nonblocked_to_blocked(structured_fs, source_block_fs, 4, 0, 7000._c_double))
  FCTEST_CHECK(check_nonblocked_to_blocked(structured_fs, source_block_fs, 4, 3, 8000._c_double))
END_TEST

! -----------------------------------------------------------------------------

TEST(test_relayout_fieldsets)
implicit none
  type(atlas_StructuredGrid) :: grid
  type(atlas_functionspace_StructuredColumns) :: structured_fs
  type(atlas_functionspace_BlockStructuredColumns) :: source_block_fs
  type(atlas_functionspace_BlockStructuredColumns) :: target_block_fs
  type(atlas_Field) :: structured_a, structured_b
  type(atlas_Field) :: blocked_a, blocked_b
  type(atlas_Field) :: nonblocked_target_a, nonblocked_target_b
  type(atlas_Field) :: blocked_target_a, blocked_target_b
  type(atlas_Field) :: blocked_relayout_a, blocked_relayout_b
  type(atlas_FieldSet) :: structured, blocked, nonblocked_target, blocked_target
  type(atlas_FieldSet) :: source_blocked, target_blocked

  grid = atlas_StructuredGrid("O8")
  structured_fs = atlas_functionspace_StructuredColumns(grid, halo=0)
  source_block_fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=5)
  target_block_fs = atlas_functionspace_BlockStructuredColumns(grid, halo=0, nproma=7)

  structured_a = structured_fs%create_field(name="structured_a", kind=atlas_real(c_double), levels=3)
  structured_b = structured_fs%create_field(name="structured_b", kind=atlas_real(c_double), levels=4, variables=2)
  blocked_a = source_block_fs%create_field(name="blocked_a", kind=atlas_real(c_double), levels=3)
  blocked_b = source_block_fs%create_field(name="blocked_b", kind=atlas_real(c_double), levels=4, variables=2)
  nonblocked_target_a = structured_fs%create_field(name="nonblocked_target_a", kind=atlas_real(c_double), levels=3)
  nonblocked_target_b = structured_fs%create_field(name="nonblocked_target_b", kind=atlas_real(c_double), levels=4, variables=2)
  blocked_target_a = source_block_fs%create_field(name="blocked_target_a", kind=atlas_real(c_double), levels=3)
  blocked_target_b = source_block_fs%create_field(name="blocked_target_b", kind=atlas_real(c_double), levels=4, variables=2)
  blocked_relayout_a = target_block_fs%create_field(name="blocked_relayout_a", kind=atlas_real(c_double), levels=3)
  blocked_relayout_b = target_block_fs%create_field(name="blocked_relayout_b", kind=atlas_real(c_double), levels=4, variables=2)

  call fill_structured_field(structured_a, 9000._c_double)
  call fill_structured_field(structured_b, 10000._c_double)
  call fill_blocked_field(source_block_fs, blocked_a, 11000._c_double)
  call fill_blocked_field(source_block_fs, blocked_b, 12000._c_double)

  structured = atlas_FieldSet()
  call structured%add(structured_a)
  call structured%add(structured_b)
  blocked = atlas_FieldSet()
  call blocked%add(blocked_a)
  call blocked%add(blocked_b)
  nonblocked_target = atlas_FieldSet()
  call nonblocked_target%add(nonblocked_target_a)
  call nonblocked_target%add(nonblocked_target_b)
  blocked_target = atlas_FieldSet()
  call blocked_target%add(blocked_target_a)
  call blocked_target%add(blocked_target_b)
  source_blocked = atlas_FieldSet()
  call source_blocked%add(blocked_a)
  call source_blocked%add(blocked_b)
  target_blocked = atlas_FieldSet()
  call target_blocked%add(blocked_relayout_a)
  call target_blocked%add(blocked_relayout_b)

  call copy_blocked_to_nonblocked(blocked, nonblocked_target)
  call copy_nonblocked_to_blocked(structured, blocked_target)
  call copy_blocked_to_blocked(source_blocked, target_blocked)

  FCTEST_CHECK_EQUAL(source_block_fs%checksum(blocked), structured_fs%checksum(nonblocked_target))
  FCTEST_CHECK_EQUAL(structured_fs%checksum(structured), source_block_fs%checksum(blocked_target))
  FCTEST_CHECK_EQUAL(source_block_fs%checksum(source_blocked), target_block_fs%checksum(target_blocked))
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
