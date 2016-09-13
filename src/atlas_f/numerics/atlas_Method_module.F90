
module atlas_Method_module

use fckit_refcounted_module, only : fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_Method

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_Method

! Purpose :
! -------
!   *Method* :
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
  procedure, public :: name => atlas_Method__name
  procedure, public :: delete => atlas_Method__delete
  procedure, public :: copy => atlas_Method__copy
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_Method__final
#endif

END TYPE atlas_Method

interface atlas_Method
  module procedure atlas_Method__cptr
end interface

!========================================================
contains
!========================================================

function atlas_Method__cptr(cptr) result(Method)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Method) :: Method
  type(c_ptr), intent(in) :: cptr
  call Method%reset_c_ptr( cptr )
end function

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_Method__final(this)
  type(atlas_Method), intent(inout) :: this
  call this%final()
end subroutine
#endif

subroutine atlas_Method__delete(this)
  use atlas_Method_c_binding
  class(atlas_Method), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Method__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_Method__delete


subroutine atlas_Method__copy(this,obj_in)
  class(atlas_Method), intent(inout) :: this
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine


function atlas_Method__name(this) result(name)
  use atlas_Method_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Method), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__Method__name(this%c_ptr())
  name = c_ptr_to_string(name_c_str)
end function

end module atlas_Method_module

