! (C) Copyright 2025 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

!> @file atlas_relayout_module.F90
!> @brief Fortran interfaces for copying Atlas fields between blocked and nonblocked layouts.
!>
!> Blocked fields are expected to come from atlas_functionspace_BlockStructuredColumns and have
!> Fortran layout `(nproma,nblk)`, `(nproma,nlev,nblk)`, or `(nproma,nvar,nlev,nblk)`.
!> Nonblocked fields are expected to come from atlas_functionspace_StructuredColumns and have
!> Fortran layout `(npoint)`, `(npoint,nlev)`, or `(npoint,nlev,nvar)`.
module atlas_Relayout_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int

implicit none

public :: copy_blocked_to_blocked
public :: copy_blocked_to_nonblocked
public :: copy_nonblocked_to_blocked

private

interface
  subroutine atlas__copy_blocked_to_blocked_field(source, target, on_device) &
      bind(C, name="atlas__copy_blocked_to_blocked_field")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__copy_blocked_to_blocked_field

  subroutine atlas__copy_blocked_to_blocked_fieldset(source, target, on_device) &
      bind(C, name="atlas__copy_blocked_to_blocked_fieldset")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__copy_blocked_to_blocked_fieldset

  subroutine atlas__copy_blocked_to_nonblocked_field(source, target, on_device) &
      bind(C, name="atlas__copy_blocked_to_nonblocked_field")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__copy_blocked_to_nonblocked_field

  subroutine atlas__copy_blocked_to_nonblocked_fieldset(source, target, on_device) &
      bind(C, name="atlas__copy_blocked_to_nonblocked_fieldset")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__copy_blocked_to_nonblocked_fieldset

  subroutine atlas__copy_nonblocked_to_blocked_field(source, target, on_device) &
      bind(C, name="atlas__copy_nonblocked_to_blocked_field")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__copy_nonblocked_to_blocked_field

  subroutine atlas__copy_nonblocked_to_blocked_fieldset(source, target, on_device) &
      bind(C, name="atlas__copy_nonblocked_to_blocked_fieldset")
    use, intrinsic :: iso_c_binding, only : c_ptr, c_int
    type(c_ptr), value :: source
    type(c_ptr), value :: target
    integer(c_int), value :: on_device
  end subroutine atlas__copy_nonblocked_to_blocked_fieldset
end interface

!> @brief Copy blocked fields to blocked fields, allowing different nproma values.
!>
!> @param source Source atlas_Field or atlas_FieldSet in blocked layout.
!> @param target Target atlas_Field or atlas_FieldSet in blocked layout.
!> @param on_device Optional execution selector.  If present and true, source and target must
!>        have valid device data and the device implementation is used.
!>
!> @pre Source and target field datatypes must match.
!> @pre Field ranks must match and must be 2, 3, or 4.
!> @pre FieldSet source and target arguments must contain the same number of fields.
interface copy_blocked_to_blocked
  module procedure copy_blocked_to_blocked_field
  module procedure copy_blocked_to_blocked_fieldset
end interface copy_blocked_to_blocked

!> @brief Copy blocked fields to nonblocked fields.
!>
!> @param source Source atlas_Field or atlas_FieldSet in blocked layout.
!> @param target Target atlas_Field or atlas_FieldSet in nonblocked layout.
!> @param on_device Optional execution selector.  If present and true, source and target must
!>        have valid device data and the device implementation is used.
!>
!> @pre Source and target field datatypes must match.
!> @pre Target field rank must be one less than source field rank.
!> @pre FieldSet source and target arguments must contain the same number of fields.
interface copy_blocked_to_nonblocked
  module procedure copy_blocked_to_nonblocked_field
  module procedure copy_blocked_to_nonblocked_fieldset
end interface copy_blocked_to_nonblocked

!> @brief Copy nonblocked fields to blocked fields.
!>
!> @param source Source atlas_Field or atlas_FieldSet in nonblocked layout.
!> @param target Target atlas_Field or atlas_FieldSet in blocked layout.
!> @param on_device Optional execution selector.  If present and true, source and target must
!>        have valid device data and the device implementation is used.
!>
!> @pre Source and target field datatypes must match.
!> @pre Target field rank must be one greater than source field rank.
!> @pre FieldSet source and target arguments must contain the same number of fields.
interface copy_nonblocked_to_blocked
  module procedure copy_nonblocked_to_blocked_field
  module procedure copy_nonblocked_to_blocked_fieldset
end interface copy_nonblocked_to_blocked

!========================================================
contains
!========================================================

subroutine copy_blocked_to_blocked_field(source, target, on_device)
  use atlas_Field_module, only : atlas_Field
  class(atlas_Field), intent(in) :: source
  class(atlas_Field), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__copy_blocked_to_blocked_field(source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine copy_blocked_to_blocked_field

subroutine copy_blocked_to_blocked_fieldset(source, target, on_device)
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_FieldSet), intent(in) :: source
  class(atlas_FieldSet), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__copy_blocked_to_blocked_fieldset(source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine copy_blocked_to_blocked_fieldset

subroutine copy_blocked_to_nonblocked_field(source, target, on_device)
  use atlas_Field_module, only : atlas_Field
  class(atlas_Field), intent(in) :: source
  class(atlas_Field), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__copy_blocked_to_nonblocked_field(source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine copy_blocked_to_nonblocked_field

subroutine copy_blocked_to_nonblocked_fieldset(source, target, on_device)
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_FieldSet), intent(in) :: source
  class(atlas_FieldSet), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__copy_blocked_to_nonblocked_fieldset(source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine copy_blocked_to_nonblocked_fieldset

subroutine copy_nonblocked_to_blocked_field(source, target, on_device)
  use atlas_Field_module, only : atlas_Field
  class(atlas_Field), intent(in) :: source
  class(atlas_Field), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__copy_nonblocked_to_blocked_field(source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine copy_nonblocked_to_blocked_field

subroutine copy_nonblocked_to_blocked_fieldset(source, target, on_device)
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_FieldSet), intent(in) :: source
  class(atlas_FieldSet), intent(inout) :: target
  logical, optional, intent(in) :: on_device
  integer(c_int) :: on_device_int

  on_device_int = 0
  if (present(on_device)) then
    if (on_device) on_device_int = 1
  endif

  call atlas__copy_nonblocked_to_blocked_fieldset(source%c_ptr(), target%c_ptr(), on_device_int)
end subroutine copy_nonblocked_to_blocked_fieldset

end module atlas_Relayout_module
