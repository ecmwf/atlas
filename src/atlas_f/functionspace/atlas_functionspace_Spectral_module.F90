
module atlas_functionspace_Spectral_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int
use atlas_c_interop, only : c_str, c_to_f_string_cptr, atlas_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Trans_module, only: atlas_Trans

implicit none

private :: c_ptr, c_int
private :: c_str, c_to_f_string_cptr, atlas_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Trans

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
    

  procedure, private :: create_glb_field_name_kind
  procedure, private :: create_glb_field_name_kind_lev
  procedure, private :: create_glb_field_kind
  procedure, private :: create_glb_field_kind_lev
  generic, public :: create_global_field => &
    & create_glb_field_name_kind, &
    & create_glb_field_name_kind_lev, &
    & create_glb_field_kind, &
    & create_glb_field_kind_lev

  procedure, public :: gather
  procedure, public :: scatter

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

function create_field_name_kind(this,name,kind) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: kind
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind(this%c_ptr(),c_str(name),kind) )
  call field%return()
end function

function create_field_name_kind_lev(this,name,kind,levels) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: levels
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind_lev(this%c_ptr(),c_str(name),kind,levels) )
  call field%return()
end function

function create_glb_field_name_kind(this,name,kind) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: kind
  field = atlas_Field( atlas__fs__Spectral__create_global_field_name_kind(this%c_ptr(),c_str(name),kind) )
  call field%return()
end function

function create_glb_field_name_kind_lev(this,name,kind,levels) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: levels
  field = atlas_Field( atlas__fs__Spectral__create_global_field_name_kind_lev(this%c_ptr(),c_str(name),kind,levels) )
  call field%return()
end function

function create_field_kind(this,kind) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  integer(c_int), intent(in) :: kind
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind(this%c_ptr(),c_str(""),kind) )
  call field%return()
end function

function create_field_kind_lev(this,kind,levels) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: levels
  field = atlas_Field( atlas__fs__Spectral__create_field_name_kind_lev(this%c_ptr(),c_str(""),kind,levels) )
  call field%return()
end function

function create_glb_field_kind(this,kind) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  integer(c_int), intent(in) :: kind
  field = atlas_Field( atlas__fs__Spectral__create_global_field_name_kind(this%c_ptr(),c_str(""),kind) )
  call field%return()
end function

function create_glb_field_kind_lev(this,kind,levels) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: levels
  field = atlas_Field( atlas__fs__Spectral__create_global_field_name_kind_lev(this%c_ptr(),c_str(""),kind,levels) )
  call field%return()
end function

subroutine gather(this,local,global)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__SpectralFunctionSpace__gather(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine scatter(this,global,local)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__SpectralFunctionSpace__scatter(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

end module atlas_functionspace_Spectral_module

