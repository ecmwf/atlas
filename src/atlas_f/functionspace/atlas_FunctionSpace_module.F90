! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_functionspace_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_field_module, only : atlas_Field
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_config_module, only : atlas_Config

implicit none

private :: fckit_owned_object
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Config

public :: atlas_FunctionSpace

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_FunctionSpace

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

  procedure, private :: create_field_args
  procedure, private :: create_field_template
  procedure, private :: deprecated_create_field_1 ! deprecated
  procedure, private :: deprecated_create_field_2 ! deprecated

  procedure, private :: halo_exchange_field
  procedure, private :: halo_exchange_fieldset

  procedure, private :: adjoint_halo_exchange_field
  procedure, private :: adjoint_halo_exchange_fieldset

  generic, public :: create_field => &
    & create_field_args, &
    & create_field_template, &
    & deprecated_create_field_1, &
    & deprecated_create_field_2

  generic, public :: halo_exchange => halo_exchange_field, halo_exchange_fieldset
  generic, public :: adjoint_halo_exchange => adjoint_halo_exchange_field, adjoint_halo_exchange_fieldset

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_FunctionSpace__final_auto
#endif

END TYPE atlas_FunctionSpace

interface atlas_FunctionSpace
  module procedure atlas_FunctionSpace__cptr
end interface

!========================================================
contains
!========================================================

function atlas_FunctionSpace__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_FunctionSpace) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_FunctionSpace__name(this) result(name)
  use atlas_functionspace_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  integer :: size
  call atlas__FunctionSpace__name(this%CPTR_PGIBUG_A, name_c_str, size )
  name = c_ptr_to_string(name_c_str)
  call c_ptr_free(name_c_str)
end function




function create_field_args(this,kind,name,levels,variables,type,alignment,global,owner,options) result(field)
  use atlas_functionspace_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  type(atlas_Field) :: field
  class(atlas_Functionspace), intent(in) :: this
  integer,          intent(in)           :: kind
  character(len=*), intent(in), optional :: name
  integer(c_int),   intent(in), optional :: levels
  integer(c_int),   intent(in), optional :: variables
  character(len=*), intent(in), optional :: type
  integer(c_int),   intent(in), optional :: alignment
  logical,          intent(in), optional :: global
  integer(c_int),   intent(in), optional :: owner
  type(atlas_Config), intent(inout), optional :: options

  type(atlas_Config) :: options_
  options_ = atlas_Config()

  if (present(options)) then
    call options_%set(options)
  endif
  call options_%set("datatype",kind)
  if( present(name)   )    call options_%set("name",name)
  if( present(owner)  )    call options_%set("owner",owner)
  if( present(global) )    call options_%set("global",global)
  if( present(levels) )    call options_%set("levels",levels)
  if( present(variables) ) call options_%set("variables",variables)
  if( present(type) )      call options_%set("type",type)
  if( present(alignment) ) call options_%set("alignment",alignment)

  field = atlas_Field( atlas__FunctionSpace__create_field( this%CPTR_PGIBUG_A, options_%CPTR_PGIBUG_B ) )

  call field%return()
  call options_%final()
end function

!------------------------------------------------------------------------------

function create_field_template(this,template,name,global,owner,options) result(field)
  use atlas_functionspace_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  type(atlas_Field) :: field
  class(atlas_Functionspace), intent(in) :: this
  type(atlas_Field), intent(in) :: template

  character(len=*), intent(in), optional :: name
  logical,          intent(in), optional :: global
  integer(c_int),   intent(in), optional :: owner
  type(atlas_Config), intent(inout), optional :: options

  type(atlas_Config) :: options_
  options_ = atlas_Config()
  if (present(options)) then
    call options_%set(options)
  endif

  if( present(name)   )    call options_%set("name",name)
  if( present(owner)  )    call options_%set("owner",owner)
  if( present(global) )    call options_%set("global",global)

  field = atlas_Field( atlas__FunctionSpace__create_field_template( &
    & this%CPTR_PGIBUG_A, template%CPTR_PGIBUG_A,options_%CPTR_PGIBUG_B) )

  call options_%final()

  call field%return()
end function

!------------------------------------------------------------------------------

subroutine halo_exchange_fieldset(this,fieldset)
  use atlas_functionspace_c_binding
  class(atlas_Functionspace), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__FunctionSpace__halo_exchange_fieldset(this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

subroutine halo_exchange_field(this,field)
  use atlas_functionspace_c_binding
  class(atlas_Functionspace), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__FunctionSpace__halo_exchange_field(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A)
end subroutine


!------------------------------------------------------------------------------

subroutine adjoint_halo_exchange_fieldset(this,fieldset)
  use atlas_functionspace_c_binding
  class(atlas_Functionspace), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__FunctionSpace__adjoint_halo_exchange_fieldset(this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

subroutine adjoint_halo_exchange_field(this,field)
  use atlas_functionspace_c_binding
  class(atlas_Functionspace), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__FunctionSpace__adjoint_halo_exchange_field(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

! Deprecated versions compatible to support IFS CY45R1

function deprecated_create_field_1(this,name,kind,levels,vars) result(field)
  use atlas_functionspace_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  type(atlas_Field) :: field
  class(atlas_Functionspace), intent(in) :: this
  character(len=*), intent(in)           :: name
  integer,          intent(in)           :: kind
  integer(c_int),   intent(in) :: levels
  integer(c_int),   intent(in) :: vars(:)

  integer(c_int) :: opt_variables

  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("datatype",kind)
  call options%set("name",name)
  call options%set("levels",levels)
  opt_variables = sum(vars)
  call options%set("variables",opt_variables)

  field = atlas_Field( atlas__FunctionSpace__create_field( this%CPTR_PGIBUG_A, options%CPTR_PGIBUG_B ) )

  call options%final()

  call field%return()
end function

function deprecated_create_field_2(this,require_name,kind,levels) result(field)
  use atlas_functionspace_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  type(atlas_Field) :: field
  class(atlas_Functionspace), intent(in) :: this
  character(len=*), intent(in) :: require_name
  integer,          intent(in) :: kind
  integer(c_int),   intent(in) :: levels

  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("datatype",kind)
  call options%set("name",require_name)
  call options%set("levels",levels)

  field = atlas_Field( atlas__FunctionSpace__create_field( this%CPTR_PGIBUG_A, options%CPTR_PGIBUG_B ) )
  call options%final()

  call field%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_FunctionSpace__final_auto(this)
  type(atlas_FunctionSpace), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_FunctionSpace__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!------------------------------------------------------------------------------


end module atlas_functionspace_module

