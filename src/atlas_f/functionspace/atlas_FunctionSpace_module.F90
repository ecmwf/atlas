
module atlas_functionspace_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_field_module, only : atlas_Field
use atlas_config_module, only : atlas_Config

implicit none

private :: fckit_owned_object
private :: atlas_Field
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

  generic, public :: create_field => &
    & create_field_args, &
    & create_field_template, &
    & deprecated_create_field_1, &
    & deprecated_create_field_2


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
  call atlas__FunctionSpace__name(this%c_ptr(), name_c_str, size )
  name = c_ptr_to_string(name_c_str)
  call c_ptr_free(name_c_str)
end function




function create_field_args(this,kind,name,levels,variables,global,owner) result(field)
  use atlas_functionspace_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  type(atlas_Field) :: field
  class(atlas_Functionspace), intent(in) :: this
  integer,          intent(in)           :: kind
  character(len=*), intent(in), optional :: name
  integer(c_int),   intent(in), optional :: levels
  integer(c_int),   intent(in), optional :: variables
  logical,          intent(in), optional :: global
  integer(c_int),   intent(in), optional :: owner

  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("datatype",kind)
  if( present(name)   )    call options%set("name",name)
  if( present(owner)  )    call options%set("owner",owner)
  if( present(global) )    call options%set("global",global)
  if( present(levels) )    call options%set("levels",levels)
  if( present(variables) ) call options%set("variables",variables)

  field = atlas_Field( atlas__FunctionSpace__create_field( this%c_ptr(), options%c_ptr() ) )
  
  call field%return()
  call options%final()
end function

!------------------------------------------------------------------------------

function create_field_template(this,template,name,global,owner) result(field)
  use atlas_functionspace_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  type(atlas_Field) :: field
  class(atlas_Functionspace), intent(in) :: this
  type(atlas_Field), intent(in) :: template

  character(len=*), intent(in), optional :: name
  logical,          intent(in), optional :: global
  integer(c_int),   intent(in), optional :: owner

  type(atlas_Config) :: options
  options = atlas_Config()

  if( present(name)   )    call options%set("name",name)
  if( present(owner)  )    call options%set("owner",owner)
  if( present(global) )    call options%set("global",global)

  field = atlas_Field( atlas__FunctionSpace__create_field_template( &
    & this%c_ptr(), template%c_ptr(),options%c_ptr()) )

  call options%final()

  call field%return()
end function

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

  field = atlas_Field( atlas__FunctionSpace__create_field( this%c_ptr(), options%c_ptr() ) )

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

  integer(c_int) :: opt_variables

  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("datatype",kind)
  call options%set("name",require_name)
  call options%set("levels",levels)

  field = atlas_Field( atlas__FunctionSpace__create_field( this%c_ptr(), options%c_ptr() ) )
  call options%final()

  call field%return()
end function

!------------------------------------------------------------------------------


end module atlas_functionspace_module

