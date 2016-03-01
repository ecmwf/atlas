
module atlas_functionspace_ReducedGridColumns_module

use, intrinsic :: iso_c_binding, only : c_ptr
use atlas_c_interop, only : c_str, c_to_f_string_cptr, atlas_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Grid_module, only: atlas_Grid

implicit none

private :: c_ptr
private :: c_str, c_to_f_string_cptr, atlas_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Grid

public :: atlas_functionspace_ReducedGridColumns

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_ReducedGridColumns

! Purpose :
! -------
!   *atlas_functionspace_ReducedGridColumns* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: create_field_name     => ReducedGridColumns__create_field_name
  procedure, private :: create_field_name_lev => ReducedGridColumns__create_field_name_lev
  generic, public :: create_field => &
    & create_field_name, &
    & create_field_name_lev

  procedure, private :: create_glb_field_name     => ReducedGridColumns__create_glb_field_name
  procedure, private :: create_glb_field_name_lev => ReducedGridColumns__create_glb_field_name_lev
  generic, public :: create_global_field => &
    & create_glb_field_name, &
    & create_glb_field_name_lev

  procedure, public :: gather => ReducedGridColumns__gather
  procedure, public :: scatter => ReducedGridColumns__scatter

  procedure, private :: checksum_fieldset => ReducedGridColumns__checksum_fieldset
  procedure, private :: checksum_field => ReducedGridColumns__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field


#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_ReducedGridColumns__final
#endif

END TYPE atlas_functionspace_ReducedGridColumns

interface atlas_functionspace_ReducedGridColumns
  module procedure ReducedGridColumns__cptr
  module procedure ReducedGridColumns__grid
end interface


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

function ReducedGridColumns__cptr(cptr) result(functionspace)
  type(atlas_functionspace_ReducedGridColumns) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function


#ifdef FORTRAN_SUPPORTS_FINAL
subroutine ReducedGridColumns__final(this)
  type(atlas_functionspace_ReducedGridColumns), intent(inout) :: this
  call this%final()
end subroutine
#endif

function ReducedGridColumns__grid(grid) result(functionspace)
  use atlas_functionspace_ReducedGridColumns_c_binding
  type(atlas_functionspace_ReducedGridColumns) :: functionspace
  class(atlas_Grid), intent(in) :: grid
  functionspace = ReducedGridColumns__cptr( &
    & atlas__functionspace__ReducedGridColumns__new__grid( grid%c_ptr()) )
  call functionspace%return()
end function

function ReducedGridColumns__create_field_name(this,name) result(field)
  use atlas_functionspace_ReducedGridColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridColumns) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__functionspace__ReducedGridColumns__create_field(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function ReducedGridColumns__create_field_name_lev(this,name,levels) result(field)
  use atlas_functionspace_ReducedGridColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridColumns), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__functionspace__ReducedGridColumns__create_field_lev(this%c_ptr(),c_str(name),levels) )
  call field%return()
end function

function ReducedGridColumns__create_glb_field_name(this,name) result(field)
  use atlas_functionspace_ReducedGridColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridColumns) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__functionspace__ReducedGridColumns__create_gfield(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function ReducedGridColumns__create_glb_field_name_lev(this,name,levels) result(field)
  use atlas_functionspace_ReducedGridColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridColumns), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__functionspace__ReducedGridColumns__create_gfield_lev(this%c_ptr(),c_str(name),levels) )
  call field%return()
end function

subroutine ReducedGridColumns__gather(this,local,global)
  use atlas_functionspace_ReducedGridColumns_c_binding
  class(atlas_functionspace_ReducedGridColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__functionspace__ReducedGridColumns__gather(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine ReducedGridColumns__scatter(this,global,local)
  use atlas_functionspace_ReducedGridColumns_c_binding
  class(atlas_functionspace_ReducedGridColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__functionspace__ReducedGridColumns__scatter(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

function ReducedGridColumns__checksum_fieldset(this,fieldset) result(checksum)
  use atlas_functionspace_ReducedGridColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_ReducedGridColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__fs__ReducedGridColumns__checksum_fieldset( &
    & this%c_ptr(),fieldset%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


function ReducedGridColumns__checksum_field(this,field) result(checksum)
  use atlas_functionspace_ReducedGridColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_ReducedGridColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__fs__ReducedGridColumns__checksum_field(this%c_ptr(),field%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function

end module atlas_functionspace_ReducedGridColumns_module

