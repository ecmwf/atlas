
module atlas_numerics_Method_module

use atlas_refcounted_module, only : atlas_RefCounted

implicit none

private :: atlas_RefCounted

public :: atlas_numerics_Method

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_numerics_Method

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
  procedure, public :: name => atlas_numerics_Method__name
  procedure, public :: delete => atlas_numerics_Method__delete
  procedure, public :: copy => atlas_numerics_Method__copy
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_numerics_Method__final
#endif

END TYPE atlas_numerics_Method

interface atlas_numerics_Method
  module procedure atlas_numerics_Method__cptr
end interface

!========================================================
contains
!========================================================

function atlas_numerics_Method__cptr(cptr) result(Method)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_numerics_Method) :: Method
  type(c_ptr), intent(in) :: cptr
  call Method%reset_c_ptr( cptr )
end function

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_numerics_Method__final(this)
  type(atlas_numerics_Method), intent(inout) :: this
  call this%final()
end subroutine
#endif

subroutine atlas_numerics_Method__delete(this)
  use atlas_numerics_Method_c_binding
  class(atlas_numerics_Method), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Method__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_numerics_Method__delete


subroutine atlas_numerics_Method__copy(this,obj_in)
  class(atlas_numerics_Method), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine


function atlas_numerics_Method__name(this) result(name)
  use atlas_numerics_Method_c_binding
  use atlas_c_interop, only : c_to_f_string_cptr
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_numerics_Method), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__Method__name(this%c_ptr())
  name = c_to_f_string_cptr(name_c_str)
end function

end module atlas_numerics_Method_module

