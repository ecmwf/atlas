#include "atlas/atlas_f.h"

module atlas_FieldSet_module

use fckit_owned_object_module, only: fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_FieldSet

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_FieldSet

! Purpose :
! -------
!   *FieldSet* : Object that groups Fields that go together
!       Fields can belong to several fieldsets simultaneously.
!       The actual ownership of the field lies in a FunctionSpace

! Methods :
! -------
!   add_field : The name or tag this field was created with
!   field : Return the field as a fortran array of specified shape

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: size => FieldSet__size
  procedure, public :: has_field
  procedure, private :: field_by_name
  procedure, private :: field_by_idx_int
  procedure, private :: field_by_idx_size_t
  procedure, public :: add
  generic :: field => field_by_name, field_by_idx_int, field_by_idx_size_t

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_FieldSet__final_auto
#endif
END TYPE atlas_FieldSet
!------------------------------------------------------------------------------

interface atlas_FieldSet
  module procedure atlas_FieldSet__cptr
  module procedure atlas_FieldSet__ctor
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! FieldSet routines

function atlas_FieldSet__cptr(cptr) result(fieldset)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_FieldSet) :: fieldset
  type(c_ptr), intent(in) :: cptr
  call fieldset%reset_c_ptr( cptr )
  call fieldset%return()
end function

function atlas_FieldSet__ctor(name) result(fieldset)
  use fckit_c_interop_module, only: c_str
  use atlas_fieldset_c_binding
  character(len=*), intent(in), optional :: name
  type(atlas_FieldSet) :: fieldset
  if( present(name) ) then
    fieldset = atlas_FieldSet__cptr( atlas__FieldSet__new( c_str(name) ) )
  else
    fieldset = atlas_FieldSet__cptr( atlas__FieldSet__new( c_str("") ) )
  endif
  call fieldset%return()
end function

subroutine add(this,field)
  use atlas_fieldset_c_binding
  use atlas_Field_module, only: atlas_Field
  class(atlas_FieldSet), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__FieldSet__add_field(this%c_ptr(), field%c_ptr())
end subroutine

function has_field(this,name) result(flag)
  use, intrinsic :: iso_c_binding, only: c_int
  use fckit_c_interop_module, only: c_str
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer(c_int) :: rc
  rc = atlas__FieldSet__has_field(this%c_ptr(), c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function

function FieldSet__size(this) result(nb_fields)
  use, intrinsic :: iso_c_binding, only: c_size_t
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  integer(c_size_t) :: nb_fields
  nb_fields = atlas__FieldSet__size(this%c_ptr())
end function

function field_by_name(this,name) result(field)
  use fckit_c_interop_module, only: c_str
  use atlas_fieldset_c_binding
  use atlas_Field_module, only: atlas_Field
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_name(this%c_ptr(), c_str(name) ) )
  call field%return()
end function

function field_by_idx_size_t(this,idx) result(field)
  use, intrinsic :: iso_c_binding, only: c_size_t
  use atlas_fieldset_c_binding
  use atlas_Field_module, only: atlas_Field
  class(atlas_FieldSet), intent(in) :: this
  integer(c_size_t), intent(in) :: idx
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_idx(this%c_ptr(), idx-1_c_size_t) ) ! C index
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use, intrinsic :: iso_c_binding, only: c_size_t, c_int
  use atlas_fieldset_c_binding
  use atlas_Field_module, only: atlas_Field
  class(atlas_FieldSet), intent(in) :: this
  integer(c_int), intent(in) :: idx
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_idx(this%c_ptr(), int(idx-1,c_size_t) ) ) ! C index
  call field%return()
end function

!-------------------------------------------------------------------------------

subroutine atlas_FieldSet__final_auto(this)
  type(atlas_FieldSet) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_FieldSet__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_FieldSet_module
