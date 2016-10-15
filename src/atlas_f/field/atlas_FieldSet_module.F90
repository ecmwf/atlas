
module atlas_FieldSet_module

use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_size_t
use fckit_c_interop_module, only: c_str
use fckit_refcounted_module, only: fckit_refcounted
use atlas_Field_module, only: atlas_Field

implicit none

private :: c_ptr, c_int, c_size_t
private :: c_str
private :: fckit_refcounted
private :: atlas_Field

public :: atlas_FieldSet

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_FieldSet

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
  procedure, public :: delete
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
  type(atlas_FieldSet) :: fieldset
  type(c_ptr), intent(in) :: cptr
  call fieldset%reset_c_ptr( cptr )
end function

function atlas_FieldSet__ctor(name) result(fieldset)
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

subroutine delete(this)
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__FieldSet__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine add(this,field)
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__FieldSet__add_field(this%c_ptr(), field%c_ptr())
end subroutine

function has_field(this,name) result(flag)
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
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  integer(c_size_t) :: nb_fields
  nb_fields = atlas__FieldSet__size(this%c_ptr())
end function

function field_by_name(this,name) result(field)
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_name(this%c_ptr(), c_str(name) ) )
  call field%return()
end function

function field_by_idx_size_t(this,idx) result(field)
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  integer(c_size_t), intent(in) :: idx
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_idx(this%c_ptr(), idx-1_c_size_t) ) ! C index
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use atlas_fieldset_c_binding
  class(atlas_FieldSet), intent(in) :: this
  integer(c_int), intent(in) :: idx
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_idx(this%c_ptr(), int(idx-1,c_size_t) ) ) ! C index
  call field%return()
end function

! ----------------------------------------------------------------------------------------

end module atlas_FieldSet_module
