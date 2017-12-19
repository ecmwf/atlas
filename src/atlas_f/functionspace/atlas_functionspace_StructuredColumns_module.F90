
module atlas_functionspace_StructuredColumns_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int
use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
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
  integer(c_int), pointer, public :: index(:,:)

contains

  procedure, public :: shared_ptr_cast

  procedure, public :: gather
  procedure, public :: scatter

  procedure, private :: checksum_fieldset
  procedure, private :: checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field

  procedure, private :: halo_exchange_fieldset
  procedure, private :: halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_fieldset, halo_exchange_field

  procedure :: j_begin
  procedure :: j_end
  procedure :: i_begin
  procedure :: i_end
  procedure :: j_begin_halo
  procedure :: j_end_halo
  procedure :: i_begin_halo
  procedure :: i_end_halo

  procedure :: xy
  procedure :: partition
  procedure :: global_index

  procedure, private :: set_index

END TYPE atlas_functionspace_StructuredColumns

interface atlas_functionspace_StructuredColumns
  module procedure StructuredColumns__cptr
  module procedure StructuredColumns__grid
end interface


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

function shared_ptr_cast(this) result(success)
  class(atlas_functionspace_StructuredColumns) :: this
  logical :: success
  success = this%fckit_shared_object__shared_ptr_cast()
  if( success ) then
    call this%set_index()
  endif
end function

function StructuredColumns__cptr(cptr) result(this)
  type(atlas_functionspace_StructuredColumns) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%set_index()
  call this%return()
end function

function StructuredColumns__grid(grid, halo, levels) result(this)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_functionspace_StructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  integer, optional :: halo
  integer, optional :: levels
  type(atlas_Config) :: config
  config = atlas_Config()
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__functionspace__StructuredColumns__new__grid( grid%c_ptr(), config%c_ptr() ) )
  call this%set_index()
  call config%final()
  call this%return()
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

subroutine set_index(this)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_f_pointer
  use atlas_functionspace_StructuredColumns_c_binding
  class(atlas_functionspace_StructuredColumns), intent(inout) :: this
  type(c_ptr) :: index_cptr
  integer(c_int), pointer :: index_fptr(:)
  integer(c_int) :: i_min, i_max, j_min, j_max
  integer(c_int) :: size
  integer(c_int) :: ni, nj
  call atlas__fs__StructuredColumns__index_host( this%c_ptr(), index_cptr, i_min, i_max, j_min, j_max )
  ni = i_max-i_min+1;
  nj = j_max-j_min+1;
  call c_f_pointer( index_cptr, index_fptr, (/ni*nj/) )
  this%index(i_min:i_max,j_min:j_max) => index_fptr(:)
end subroutine

function j_begin(this) result(j)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_begin(this%c_ptr())
end function

function j_end(this) result(j)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_end(this%c_ptr())
end function

function i_begin(this,j) result(i)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: i
  integer(c_int), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_begin(this%c_ptr(),j)
end function

function i_end(this,j) result(i)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: i
  integer(c_int), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_end(this%c_ptr(),j)
end function


function j_begin_halo(this) result(j)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_begin_halo(this%c_ptr())
end function

function j_end_halo(this) result(j)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_end_halo(this%c_ptr())
end function

function i_begin_halo(this,j) result(i)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: i
  integer(c_int), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_begin_halo(this%c_ptr(),j)
end function

function i_end_halo(this,j) result(i)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_functionspace_StructuredColumns_c_binding
  integer(c_int) :: i
  integer(c_int), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_end_halo(this%c_ptr(),j)
end function

function xy(this) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__StructuredColumns__xy(this%c_ptr()) )
  call field%return()
end function

function partition(this) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__StructuredColumns__partition(this%c_ptr()) )
  call field%return()
end function


function global_index(this) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__StructuredColumns__global_index(this%c_ptr()) )
  call field%return()
end function

subroutine halo_exchange_fieldset(this,fieldset)
  use atlas_functionspace_StructuredColumns_c_binding
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__fs__StructuredColumns__halo_exchange_fieldset(this%c_ptr(),fieldset%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine halo_exchange_field(this,field)
  use atlas_functionspace_StructuredColumns_c_binding
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__fs__StructuredColumns__halo_exchange_field(this%c_ptr(),field%c_ptr())
end subroutine


end module atlas_functionspace_StructuredColumns_module

