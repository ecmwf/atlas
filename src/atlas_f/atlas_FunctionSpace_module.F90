
module atlas_functionspace_module

use atlas_refcounted_module, only : atlas_RefCounted

implicit none

private :: atlas_RefCounted

public :: atlas_NextFunctionSpace

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_NextFunctionSpace

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
  procedure, public :: name => atlas_NextFunctionSpace__name
  procedure, public :: delete => atlas_NextFunctionSpace__delete
  procedure, public :: copy => atlas_NextFunctionSpace__copy
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_NextFunctionSpace__final
#endif

END TYPE atlas_NextFunctionSpace

interface atlas_NextFunctionSpace
  module procedure atlas_NextFunctionSpace__cptr
end interface

!========================================================
contains
!========================================================

function atlas_NextFunctionSpace__cptr(cptr) result(functionspace)
  use iso_c_binding, only : c_ptr
  type(atlas_NextFunctionSpace) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_NextFunctionSpace__final(this)
  type(atlas_NextFunctionSpace), intent(inout) :: this
  call this%final()
end subroutine
#endif

subroutine atlas_NextFunctionSpace__delete(this)
  use atlas_functionspace_c_binding
  class(atlas_NextFunctionSpace), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__NextFunctionSpace__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_NextFunctionSpace__delete


subroutine atlas_NextFunctionSpace__copy(this,obj_in)
  class(atlas_NextFunctionSpace), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine


function atlas_NextFunctionSpace__name(this) result(name)
  use atlas_functionspace_c_binding
  use atlas_c_interop, only : c_to_f_string_cptr
  use iso_c_binding, only : c_ptr
  class(atlas_NextFunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__NextFunctionSpace__name(this%c_ptr())
  name = c_to_f_string_cptr(name_c_str)
end function

end module atlas_functionspace_module

