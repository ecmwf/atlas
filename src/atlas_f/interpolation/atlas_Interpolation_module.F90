#include "atlas/atlas_f.h"

module atlas_Interpolation_module

use fckit_owned_object_module, only : fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_Interpolation

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Interpolation

! Purpose :
! -------
!   *Interpolation* :

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: execute_field
  procedure, private :: execute_fieldset
  generic, public :: execute => execute_field, execute_fieldset

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Interpolation__final_auto
#endif

END TYPE atlas_Interpolation

interface atlas_Interpolation
  module procedure atlas_Interpolation__cptr
  module procedure atlas_Interpolation__config_funcspace
end interface

!========================================================
contains
!========================================================

function atlas_Interpolation__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_Interpolation) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function

function atlas_Interpolation__config_funcspace(config,source,target) result(this)
  use atlas_Interpolation_c_binding
  use atlas_Config_module, only : atlas_Config
  use atlas_FunctionSpace_module, only : atlas_FunctionSpace
  type(atlas_Interpolation) :: this
  type(atlas_Config), intent(in) :: config
  class(atlas_FunctionSpace), intent(in) :: source
  class(atlas_FunctionSpace), intent(in) :: target
  this = atlas_Interpolation__cptr(atlas__interpolation__new(config%c_ptr(),source%c_ptr(),target%c_ptr()))
  call this%return()
end function

subroutine execute_field(this,source,target)
  use atlas_Interpolation_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Interpolation), intent(in) :: this
  class(atlas_Field), intent(in) :: source
  class(atlas_Field), intent(inout) :: target
  call atlas__Interpolation__execute_field(this%c_ptr(),source%c_ptr(),target%c_ptr())
end subroutine

subroutine execute_fieldset(this,source,target)
  use atlas_Interpolation_c_binding
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_Interpolation), intent(in) :: this
  class(atlas_FieldSet), intent(in) :: source
  class(atlas_FieldSet), intent(inout) :: target
  call atlas__Interpolation__execute_fieldset(this%c_ptr(),source%c_ptr(),target%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

subroutine atlas_Interpolation__final_auto(this)
  type(atlas_Interpolation) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Interpolation__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! -----------------------------------------------------------------------------

end module atlas_Interpolation_module

