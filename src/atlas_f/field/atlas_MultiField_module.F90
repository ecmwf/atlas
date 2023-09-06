! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_multifield_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_Config_module, only: atlas_Config
use atlas_field_module, only: atlas_field, array_c_to_f
use atlas_fieldset_module, only: atlas_fieldset

implicit none

public :: atlas_MultiField

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_MultiField

! Purpose :
! -------
!   *MultiField* : Object containing field data and Metadata

! Methods :
! -------
!   name : The name or tag this field was created with
!   data : Return the field as a fortran array of specified shape
!   Metadata : Return object that can contain a variety of Metadata

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*
!   29-Aug-2023 Slavko Brdar         *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public  :: MultiField__fieldset
  procedure, public  :: MultiField__size
  generic :: fieldset => MultiField__fieldset
  generic :: size => MultiField__size

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_MultiField__final_auto
#endif

END TYPE

interface atlas_MultiField
  module procedure atlas_MultiField__cptr
  module procedure atlas_MultiField__create
  module procedure atlas_MultiField__create_names
end interface

private :: fckit_owned_object
private :: atlas_Config

!========================================================
contains
!========================================================

!-------------------------------------------------------------------------------

function atlas_MultiField__cptr(cptr) result(field)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_MultiField) :: field
  type(c_ptr), intent(in) :: cptr
  call field%reset_c_ptr( cptr )
  call field%return()
end function

!-------------------------------------------------------------------------------

function atlas_MultiField__create(params) result(field)
  use atlas_multifield_c_binding
  type(atlas_MultiField) :: field
  class(atlas_Config), intent(in) :: params
  field = atlas_MultiField__cptr( atlas__MultiField__create(params%CPTR_PGIBUG_B) )
  call field%return()
end function

!-------------------------------------------------------------------------------

function atlas_MultiField__create_names(kind, shape, field_names) result(field)
  use, intrinsic :: iso_c_binding, only : c_char, c_int, c_int32_t, c_size_t
  use atlas_multifield_c_binding
  type(atlas_MultiField) :: field
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: shape(:)
  character(*), intent(in) :: field_names(:)
  character(kind=c_char,len=:), allocatable :: flat_field_names
  integer(c_size_t) :: length
  integer(c_int32_t) :: ii
  integer(c_int32_t), allocatable :: field_name_sizes(:)

  if (size(field_names) == 0 .or. size(shape) < 3) then
    print *, "atlas_MultiField__create_names: must have at least one field name, and the size of shape", &
        & " is minimum 3, e.g. [nproma,-1,nblk]"
    stop -1
  end if

  length = len(field_names(1))
  allocate(field_name_sizes(size(field_names)))
  field_name_sizes = len(field_names(:))

  if (any(field_name_sizes /= length)) then
    print *, "atlas_MultiField__create_names: field_names have to have same length in characters"
    stop -1
  end if

  allocate(character(len=length*size(field_names) ) :: flat_field_names)
  do ii = 1, size(field_names)
    flat_field_names((ii-1)*length+1:ii*length) = field_names(ii)
  enddo

  field = atlas_MultiField__cptr( atlas__MultiField__create_shape(kind, size(shape), shape, &
        & flat_field_names, length, size(field_names,kind=c_size_t)) )
  call field%return()
end function

!-------------------------------------------------------------------------------

function MultiField__size(this) result(size)
  use atlas_multifield_c_binding
  class(atlas_MultiField), intent(in) :: this
  integer :: size
  size = atlas__MultiField__size(this%CPTR_PGIBUG_B)
end function

!-------------------------------------------------------------------------------
function MultiField__fieldset(this) result(fset)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_multifield_c_binding
  class(atlas_MultiField), intent(in) :: this
  type(c_ptr) :: fset_cptr
  type(atlas_FieldSet) :: fset
  fset_cptr = atlas__MultiField__fieldset(this%CPTR_PGIBUG_B)
  fset = atlas_FieldSet( fset_cptr )
  call fset%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_MultiField__final_auto(this)
  type(atlas_MultiField), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_MultiField__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!-------------------------------------------------------------------------------

end module atlas_multifield_module

