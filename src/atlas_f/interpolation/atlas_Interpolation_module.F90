! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

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
  procedure, private :: execute_adjoint_field
  procedure, private :: execute_adjoint_fieldset
  generic, public :: execute => execute_field, execute_fieldset
  generic, public :: execute_adjoint => execute_adjoint_field, execute_adjoint_fieldset

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Interpolation__final_auto
#endif

END TYPE atlas_Interpolation

interface atlas_Interpolation
  module procedure atlas_Interpolation__cptr
  module procedure atlas_Interpolation__config_funcspace
  module procedure atlas_Interpolation__config_funcspace_field
  module procedure atlas_Interpolation__config_funcspace_fieldset
end interface

!========================================================
contains
!========================================================

function atlas_Interpolation__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_Interpolation) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_Interpolation__config_funcspace(config,source,target) result(this)
  use atlas_Interpolation_c_binding
  use atlas_Config_module, only : atlas_Config
  use atlas_FunctionSpace_module, only : atlas_FunctionSpace
  type(atlas_Interpolation) :: this
  type(atlas_Config), intent(in) :: config
  class(atlas_FunctionSpace), intent(in) :: source
  class(atlas_FunctionSpace), intent(in) :: target
  this = atlas_Interpolation__cptr(atlas__interpolation__new(config%CPTR_PGIBUG_B, &
      source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A))
  call this%return()
end function

function atlas_Interpolation__config_funcspace_field(config,source,target) result(this)
  use atlas_Interpolation_c_binding
  use atlas_Config_module, only : atlas_Config
  use atlas_FunctionSpace_module, only : atlas_FunctionSpace
  use atlas_Field_module, only : atlas_Field
  type(atlas_Interpolation) :: this
  type(atlas_Config), intent(in) :: config
  class(atlas_FunctionSpace), intent(in) :: source
  class(atlas_Field), intent(in) :: target
  this = atlas_Interpolation__cptr(atlas__interpolation__new_tgt_field(config%CPTR_PGIBUG_B, &
      source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A))
  call this%return()
end function

function atlas_Interpolation__config_funcspace_fieldset(config,source,target) result(this)
  use atlas_Interpolation_c_binding
  use atlas_Config_module, only : atlas_Config
  use atlas_FunctionSpace_module, only : atlas_FunctionSpace
  use atlas_FieldSet_module, only : atlas_FieldSet
  type(atlas_Interpolation) :: this
  type(atlas_Config), intent(in) :: config
  class(atlas_FunctionSpace), intent(in) :: source
  class(atlas_FieldSet), intent(in) :: target
  this = atlas_Interpolation__cptr(atlas__interpolation__new_tgt_fieldset(config%CPTR_PGIBUG_B, &
      source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A))
  call this%return()
end function

subroutine execute_field(this,source,target)
  use atlas_Interpolation_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Interpolation), intent(in) :: this
  class(atlas_Field), intent(in) :: source
  class(atlas_Field), intent(inout) :: target
  call atlas__Interpolation__execute_field(this%CPTR_PGIBUG_A,source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A)
end subroutine

subroutine execute_fieldset(this,source,target)
  use atlas_Interpolation_c_binding
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_Interpolation), intent(in) :: this
  class(atlas_FieldSet), intent(in) :: source
  class(atlas_FieldSet), intent(inout) :: target
  call atlas__Interpolation__execute_fieldset(this%CPTR_PGIBUG_A,source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A)
end subroutine

subroutine execute_adjoint_field(this,source,target)
  use atlas_Interpolation_c_binding
  use atlas_Field_module, only : atlas_Field
  class(atlas_Interpolation), intent(in) :: this
  class(atlas_Field), intent(inout) :: source
  class(atlas_Field), intent(in) :: target
  call atlas__Interpolation__execute_adjoint_field(this%CPTR_PGIBUG_A,source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A)
end subroutine

subroutine execute_adjoint_fieldset(this,source,target)
  use atlas_Interpolation_c_binding
  use atlas_FieldSet_module, only : atlas_FieldSet
  class(atlas_Interpolation), intent(in) :: this
  class(atlas_FieldSet), intent(inout) :: source
  class(atlas_FieldSet), intent(in) :: target
  call atlas__Interpolation__execute_adjoint_fieldset(this%CPTR_PGIBUG_A,source%CPTR_PGIBUG_A,target%CPTR_PGIBUG_A)
end subroutine


!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Interpolation__final_auto(this)
  type(atlas_Interpolation), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Interpolation__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! -----------------------------------------------------------------------------

end module atlas_Interpolation_module

