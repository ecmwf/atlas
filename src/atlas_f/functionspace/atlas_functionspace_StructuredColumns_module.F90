
module atlas_functionspace_StructuredColumns_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int
use fckit_c_interop, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Grid_module, only: atlas_Grid
use atlas_Config_module, only: atlas_Config

implicit none

private :: c_ptr, c_int
private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Grid
private :: atlas_Config

public :: atlas_functionspace_StructuredColumns

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_StructuredColumns

! Purpose :
! -------
!   *atlas_functionspace_StructuredColumns* : Interpretes spectral fields

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

  procedure, public :: gather
  procedure, public :: scatter

  procedure, private :: checksum_fieldset
  procedure, private :: checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field


#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_StructuredColumns__final
#endif

END TYPE atlas_functionspace_StructuredColumns

interface atlas_functionspace_StructuredColumns
  module procedure StructuredColumns__cptr
  module procedure StructuredColumns__grid
end interface


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

function StructuredColumns__cptr(cptr) result(functionspace)
  type(atlas_functionspace_StructuredColumns) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function


#ifdef FORTRAN_SUPPORTS_FINAL
subroutine StructuredColumns__final(this)
  type(atlas_functionspace_StructuredColumns), intent(inout) :: this
  call this%final()
end subroutine
#endif

function StructuredColumns__grid(grid) result(functionspace)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_functionspace_StructuredColumns) :: functionspace
  class(atlas_Grid), intent(in) :: grid
  functionspace = StructuredColumns__cptr( &
    & atlas__functionspace__StructuredColumns__new__grid( grid%c_ptr()) )
  call functionspace%return()
end function

function create_field_name_kind(this,name,kind,global,owner) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns) :: this
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
  field = atlas_Field( atlas__fs__StructuredColumns__create_field_name_kind( &
      & this%c_ptr(),c_str(name),kind,options%c_ptr()) )
  call field%return()
  call options%final()
end function

function create_field_name_kind_lev(this,name,kind,levels,global,owner) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
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
  field = atlas_Field( atlas__fs__StructuredColumns__create_field_name_kind_lev( &
      & this%c_ptr(),c_str(name),kind,levels,options%c_ptr()) )
  call field%return()
  call options%final()
end function

function create_field_kind(this,kind,global,owner) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns) :: this
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
  field = atlas_Field( atlas__fs__StructuredColumns__create_field_name_kind( &
      & this%c_ptr(),c_str(""),kind,options%c_ptr()) )
  call field%return()
  call options%final()
end function

function create_field_kind_lev(this,kind,levels,global,owner) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
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
  field = atlas_Field( atlas__fs__StructuredColumns__create_field_name_kind_lev( &
      & this%c_ptr(),c_str(""),kind,levels,options%c_ptr()) )
  call field%return()
  call options%final()
end function

subroutine gather(this,local,global)
  use atlas_functionspace_StructuredColumns_c_binding
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__functionspace__StructuredColumns__gather(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine scatter(this,global,local)
  use atlas_functionspace_StructuredColumns_c_binding
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__functionspace__StructuredColumns__scatter(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

function checksum_fieldset(this,fieldset) result(checksum)
  use atlas_functionspace_StructuredColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__fs__StructuredColumns__checksum_fieldset( &
    & this%c_ptr(),fieldset%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function


function checksum_field(this,field) result(checksum)
  use atlas_functionspace_StructuredColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__fs__StructuredColumns__checksum_field( &
      & this%c_ptr(),field%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

end module atlas_functionspace_StructuredColumns_module

