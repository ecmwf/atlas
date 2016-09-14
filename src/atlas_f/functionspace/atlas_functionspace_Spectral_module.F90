
module atlas_functionspace_Spectral_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int
use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Trans_module, only: atlas_Trans
use atlas_Config_module, only: atlas_Config

implicit none

private :: c_ptr, c_int
private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Trans
private :: atlas_Config

public :: atlas_functionspace_Spectral

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_Spectral

! Purpose :
! -------
!   *atlas_functionspace_Spectral* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: create_field_name_kind
  procedure, private :: create_field_name_kind_lev
  procedure, private :: create_field_kind
  procedure, private :: create_field_kind_lev
  generic, public :: create_field => &
    & create_field_name_kind, &
    & create_field_name_kind_lev, &
    & create_field_kind, &
    & create_field_kind_lev

  procedure, private :: gather_field
  procedure, private :: scatter_field
  procedure, private :: gather_fieldset
  procedure, private :: scatter_fieldset

  generic, public :: gather  =>  gather_field,  gather_fieldset
  generic, public :: scatter => scatter_field, scatter_fieldset


#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_Spectral__final
#endif

END TYPE atlas_functionspace_Spectral

interface atlas_functionspace_Spectral
  module procedure atlas_functionspace_Spectral__cptr
  module procedure atlas_functionspace_Spectral__truncation
  module procedure atlas_functionspace_Spectral__trans
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

function atlas_functionspace_Spectral__cptr(cptr) result(functionspace)
  type(atlas_functionspace_Spectral) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function


#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_functionspace_Spectral__final(this)
  type(atlas_functionspace_Spectral), intent(inout) :: this
  call this%final()
end subroutine
#endif

function atlas_functionspace_Spectral__truncation(truncation) result(functionspace)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: functionspace
  integer(c_int), intent(in) :: truncation
  functionspace = atlas_functionspace_Spectral__cptr( &
    & atlas__SpectralFunctionSpace__new__truncation(truncation) )
  call functionspace%return()
end function

function atlas_functionspace_Spectral__trans(trans) result(functionspace)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: functionspace
  type(atlas_Trans), intent(in) :: trans
  functionspace = atlas_functionspace_Spectral__cptr( &
    & atlas__SpectralFunctionSpace__new__trans(trans%c_ptr()) )
  call functionspace%return()
end function

function create_field_name_kind(this,name,kind,global,owner) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: kind
  logical, optional, intent(in) :: global
  integer(c_int), optional, intent(in) :: owner
  logical :: opt_global
  integer(c_int) :: opt_owner
  type(atlas_Config) :: options
  opt_owner = 0
  if( present(owner) ) opt_owner = owner
  opt_global = .false.
  if( present(global) ) opt_global = global
  options = atlas_Config()
  call options%set("global",opt_global)
  call options%set("owner",opt_owner)
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind(this%c_ptr(),c_str(name),kind,options%c_ptr()) )
  call field%return()
  call options%final()
end function

function create_field_name_kind_lev(this,name,kind,levels,global,owner) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: levels
  logical, optional, intent(in) :: global
  integer(c_int), optional, intent(in) :: owner
  logical :: opt_global
  integer(c_int) :: opt_owner
  type(atlas_Config) :: options
  opt_owner = 0
  if( present(owner) ) opt_owner = owner
  opt_global = .false.
  if( present(global) ) opt_global = global
  options = atlas_Config()
  call options%set("global",opt_global)
  call options%set("owner",opt_owner)
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind_lev(this%c_ptr(),c_str(name),kind,levels,options%c_ptr()) )
  call field%return()
  call options%final()
end function

function create_field_kind(this,kind,global,owner) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  integer(c_int), intent(in) :: kind
  logical, optional, intent(in) :: global
  integer(c_int), optional, intent(in) :: owner
  logical :: opt_global
  integer(c_int) :: opt_owner
  type(atlas_Config) :: options
  opt_owner = 0
  if( present(owner) ) opt_owner = owner
  opt_global = .false.
  if( present(global) ) opt_global = global
  options = atlas_Config()
  call options%set("global",opt_global)
  call options%set("owner",opt_owner)
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind(this%c_ptr(),c_str(""),kind,options%c_ptr()) )
  call field%return()
  call options%final()
end function

function create_field_kind_lev(this,kind,levels,global,owner) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: levels
  logical, optional, intent(in) :: global
  integer(c_int), optional, intent(in) :: owner
  logical :: opt_global
  integer(c_int) :: opt_owner
  type(atlas_Config) :: options
  opt_owner = 0
  if( present(owner) ) opt_owner = owner
  opt_global = .false.
  if( present(global) ) opt_global = global
  options = atlas_Config()
  call options%set("global",opt_global)
  call options%set("owner",opt_owner)
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind_lev(this%c_ptr(),c_str(""),kind,levels,options%c_ptr()) )
  call field%return()
  call options%final()
end function

subroutine gather_field(this,local,global)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__SpectralFunctionSpace__gather(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine scatter_field(this,global,local)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__SpectralFunctionSpace__scatter(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

subroutine gather_fieldset(this,local,global)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__SpectralFunctionSpace__gather_fieldset(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine scatter_fieldset(this,global,local)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__SpectralFunctionSpace__scatter_fieldset(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

end module atlas_functionspace_Spectral_module

