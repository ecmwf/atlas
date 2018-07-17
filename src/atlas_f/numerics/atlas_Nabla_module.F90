#include "atlas/atlas_f.h"

module atlas_Nabla_module

use fckit_owned_object_module, only : fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_Nabla

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Nabla

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
  procedure, public :: gradient => atlas_Nabla__gradient
  procedure, public :: divergence => atlas_Nabla__divergence
  procedure, public :: curl => atlas_Nabla__curl
  procedure, public :: laplacian => atlas_Nabla__laplacian

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Nabla__final_auto
#endif
END TYPE atlas_Nabla

interface atlas_Nabla
  module procedure atlas_Nabla__cptr
  module procedure atlas_Nabla__method_config
end interface

!========================================================
contains
!========================================================

function atlas_Nabla__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_Nabla) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_Nabla__method_config(method,config) result(this)
  use atlas_Nabla_c_binding
  use atlas_Method_module, only : atlas_Method
  use atlas_Config_module, only : atlas_Config
  type(atlas_Nabla) :: this
  class(atlas_Method), intent(in) :: method
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call this%reset_c_ptr( atlas__Nabla__create(method%c_ptr(),config%c_ptr()) )
  else
    opt_config = atlas_Config()
    call this%reset_c_ptr( atlas__Nabla__create(method%c_ptr(),opt_config%c_ptr()) )
    call opt_config%final()
  endif
  call this%return()
end function

subroutine atlas_Nabla__gradient(this,scalar,grad)
  use atlas_Nabla_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: scalar
  class(atlas_Field), intent(inout) :: grad
  call atlas__Nabla__gradient(this%c_ptr(),scalar%c_ptr(),grad%c_ptr())
end subroutine


subroutine atlas_Nabla__divergence(this,vector,div)
  use atlas_Nabla_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: vector
  class(atlas_Field), intent(inout) :: div
  call atlas__Nabla__divergence(this%c_ptr(),vector%c_ptr(),div%c_ptr())
end subroutine

subroutine atlas_Nabla__curl(this,vector,curl)
  use atlas_Nabla_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: vector
  class(atlas_Field), intent(inout) :: curl
  call atlas__Nabla__curl(this%c_ptr(),vector%c_ptr(),curl%c_ptr())
end subroutine

subroutine atlas_Nabla__laplacian(this,scalar,lapl)
  use atlas_Nabla_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: scalar
  class(atlas_Field), intent(inout) :: lapl
  call atlas__Nabla__laplacian(this%c_ptr(),scalar%c_ptr(),lapl%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_Nabla__final_auto(this)
  type(atlas_Nabla), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Nabla__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! -----------------------------------------------------------------------------

end module atlas_Nabla_module

