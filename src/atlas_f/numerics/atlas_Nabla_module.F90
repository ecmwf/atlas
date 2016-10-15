
module atlas_Nabla_module

use fckit_refcounted_module, only : fckit_refcounted
use atlas_Method_module, only : atlas_Method
use atlas_Field_module, only : atlas_Field
use atlas_Config_module, only : atlas_Config

implicit none

private :: fckit_refcounted
private :: atlas_Method
private :: atlas_Field
private :: atlas_Config

public :: atlas_Nabla

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_Nabla

! Purpose :
! -------
!   *Nabla* :

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_Nabla__delete
  procedure, public :: copy => atlas_Nabla__copy
  procedure, public :: gradient => atlas_Nabla__gradient
  procedure, public :: divergence => atlas_Nabla__divergence
  procedure, public :: curl => atlas_Nabla__curl
  procedure, public :: laplacian => atlas_Nabla__laplacian

END TYPE atlas_Nabla

interface atlas_Nabla
  module procedure atlas_Nabla__cptr
  module procedure atlas_Nabla__method_config
end interface

!========================================================
contains
!========================================================

function atlas_Nabla__cptr(cptr) result(nabla)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_Nabla) :: nabla
  type(c_ptr), intent(in) :: cptr
  call nabla%reset_c_ptr( cptr )
end function

function atlas_Nabla__method_config(method,config) result(nabla)
  use atlas_Nabla_c_binding
  type(atlas_Nabla) :: nabla
  class(atlas_Method), intent(in) :: method
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    nabla = atlas_Nabla__cptr(atlas__Nabla__create(method%c_ptr(),config%c_ptr()))
  else
    opt_config = atlas_Config()
    nabla = atlas_Nabla__cptr(atlas__Nabla__create(method%c_ptr(),opt_config%c_ptr()))
    call opt_config%final()
  endif
  call nabla%return()
end function

subroutine atlas_Nabla__delete(this)
  use atlas_Nabla_c_binding
  class(atlas_Nabla), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Nabla__delete(this%c_ptr())
  endif
  call this%reset_c_ptr()
end subroutine atlas_Nabla__delete

subroutine atlas_Nabla__copy(this,obj_in)
  class(atlas_Nabla), intent(inout) :: this
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine


subroutine atlas_Nabla__gradient(this,scalar,grad)
  use atlas_Nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: scalar
  class(atlas_Field), intent(inout) :: grad
  call atlas__Nabla__gradient(this%c_ptr(),scalar%c_ptr(),grad%c_ptr())
end subroutine


subroutine atlas_Nabla__divergence(this,vector,div)
  use atlas_Nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: vector
  class(atlas_Field), intent(inout) :: div
  call atlas__Nabla__divergence(this%c_ptr(),vector%c_ptr(),div%c_ptr())
end subroutine

subroutine atlas_Nabla__curl(this,vector,curl)
  use atlas_Nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: vector
  class(atlas_Field), intent(inout) :: curl
  call atlas__Nabla__curl(this%c_ptr(),vector%c_ptr(),curl%c_ptr())
end subroutine

subroutine atlas_Nabla__laplacian(this,scalar,lapl)
  use atlas_Nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: scalar
  class(atlas_Field), intent(inout) :: lapl
  call atlas__Nabla__laplacian(this%c_ptr(),scalar%c_ptr(),lapl%c_ptr())
end subroutine

! -----------------------------------------------------------------------------

end module atlas_Nabla_module

