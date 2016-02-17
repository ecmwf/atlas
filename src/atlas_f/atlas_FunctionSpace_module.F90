
module atlas_functionspace_module

use atlas_refcounted_module, only : atlas_RefCounted

implicit none

private :: atlas_RefCounted

public :: atlas_FunctionSpace

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_FunctionSpace

! Purpose :
! -------
!   *FunctionSpace* :
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done

! Methods :
! -------
!   name : The name or tag this function space was created with

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: name => atlas_FunctionSpace__name
  procedure, public :: delete => atlas_FunctionSpace__delete
  procedure, public :: copy => atlas_FunctionSpace__copy
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_FunctionSpace__final
#endif

END TYPE atlas_FunctionSpace

interface atlas_FunctionSpace
  module procedure atlas_FunctionSpace__cptr
end interface

!========================================================
contains
!========================================================

function atlas_FunctionSpace__cptr(cptr) result(functionspace)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_FunctionSpace) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_FunctionSpace__final(this)
  type(atlas_FunctionSpace), intent(inout) :: this
  call this%final()
end subroutine
#endif

subroutine atlas_FunctionSpace__delete(this)
  use atlas_functionspace_c_binding
  class(atlas_FunctionSpace), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__FunctionSpace__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_FunctionSpace__delete


subroutine atlas_FunctionSpace__copy(this,obj_in)
  class(atlas_FunctionSpace), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine


function atlas_FunctionSpace__name(this) result(name)
  use atlas_functionspace_c_binding
  use atlas_c_interop, only : c_to_f_string_cptr
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__FunctionSpace__name(this%c_ptr())
  name = c_to_f_string_cptr(name_c_str)
end function

end module atlas_functionspace_module

