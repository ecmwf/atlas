
module atlas_functionspace_module

use fckit_refcounted_module, only : fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_FunctionSpace

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_FunctionSpace

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

subroutine atlas_FunctionSpace__delete(this)
  use atlas_functionspace_c_binding
  class(atlas_FunctionSpace), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__FunctionSpace__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_FunctionSpace__delete

function atlas_FunctionSpace__name(this) result(name)
  use atlas_functionspace_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__FunctionSpace__name(this%c_ptr())
  name = c_ptr_to_string(name_c_str)
end function

end module atlas_functionspace_module

