
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

  procedure, private :: create_field_name_kind_lev
  generic, public :: create_field => &
    & create_field_name_kind_lev

  procedure, private :: gather_field
  procedure, private :: scatter_field
  procedure, private :: gather_fieldset
  procedure, private :: scatter_fieldset

  generic, public :: gather  =>  gather_field,  gather_fieldset
  generic, public :: scatter => scatter_field, scatter_fieldset

  procedure, private :: norm_scalar
  procedure, private :: norm_array
  generic, public :: norm => norm_scalar, norm_array

#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_Spectral__final
#endif

END TYPE atlas_functionspace_Spectral

interface atlas_functionspace_Spectral
  module procedure atlas_functionspace_Spectral__cptr
  module procedure atlas_functionspace_Spectral__config
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

function atlas_functionspace_Spectral__config(truncation,levels) result(this)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: this
  integer(c_int), intent(in)           :: truncation
  integer(c_int), intent(in), optional :: levels
  
  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("truncation",truncation)
  if( present(levels) ) call options%set("levels",levels)
  
  this = atlas_functionspace_Spectral__cptr( &
    & atlas__SpectralFunctionSpace__new__config(options%c_ptr()) )
  call options%final()
    
  call this%return()
end function

function atlas_functionspace_Spectral__trans(trans,levels) result(this)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: this
  type(atlas_Trans), intent(in) :: trans
  integer(c_int), intent(in), optional :: levels
  
  type(atlas_Config) :: options
  options = atlas_Config()

  if( present(levels) ) call options%set("levels",levels)

  this = atlas_functionspace_Spectral__cptr( &
    & atlas__SpectralFunctionSpace__new__trans(trans%c_ptr(), options%c_ptr() ) )
  call options%final()
  
  call this%return()
end function


function create_field_name_kind_lev(this,kind,name,levels,global,owner) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  integer(c_int),   intent(in) :: kind
  character(len=*), intent(in), optional :: name
  integer(c_int),   intent(in), optional :: levels
  logical,          intent(in), optional :: global
  integer(c_int),   intent(in), optional :: owner

  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("datatype",kind)
  if( present(name) )   call options%set("name",name)
  if( present(levels) ) call options%set("levels",levels)
  if( present(owner) )  call options%set("owner",owner)
  if( present(global) ) call options%set("global",global)
  field = atlas_Field( atlas__fs__Spectral__create_field(this%c_ptr(),options%c_ptr()) )
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

subroutine norm_scalar(this,field,norm,rank)
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  real(c_double), intent(out) :: norm
  integer(c_int), optional :: rank
  integer :: opt_rank
  real(c_double) :: norm_array(1)
  opt_rank = 0
  if( present(rank) ) opt_rank = rank
  call atlas__SpectralFunctionSpace__norm(this%c_ptr(),field%c_ptr(),norm_array,opt_rank)
  norm = norm_array(1)
end subroutine

subroutine norm_array(this,field,norm,rank)
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  real(c_double), intent(inout) :: norm(:)
  integer(c_int), optional :: rank
  integer :: opt_rank
  opt_rank = 0
  if( present(rank) ) opt_rank = rank
  call atlas__SpectralFunctionSpace__norm(this%c_ptr(),field%c_ptr(),norm,opt_rank)
end subroutine

end module atlas_functionspace_Spectral_module

