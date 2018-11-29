#include "atlas/atlas_f.h"

module atlas_functionspace_StructuredColumns_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_double
use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Grid_module, only: atlas_Grid
use atlas_Config_module, only: atlas_Config
use atlas_kinds_module, only : ATLAS_KIND_IDX
use fckit_owned_object_module, only : fckit_owned_object
use atlas_GridDistribution_module, only : atlas_GridDistribution
use atlas_Vertical_module, only : atlas_Vertical

implicit none

private :: c_ptr, c_double
private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Grid
private :: atlas_GridDistribution
private :: atlas_Vertical
private :: atlas_Config
private :: fckit_owned_object

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
  integer(ATLAS_KIND_IDX), pointer, public :: index(:,:)

contains

  procedure, public :: assignment_operator_hook

  procedure, public :: gather
  procedure, public :: scatter

  procedure, private :: checksum_fieldset
  procedure, private :: checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field

  procedure :: j_begin
  procedure :: j_end
  procedure :: i_begin
  procedure :: i_end
  procedure :: j_begin_halo
  procedure :: j_end_halo
  procedure :: i_begin_halo
  procedure :: i_end_halo

  procedure :: size => get_size
  procedure :: size_owned => get_size_owned

  procedure :: xy
    !! Return xy coordinate field
  procedure :: partition
    !! Return partition field
  procedure :: global_index
    !! Return global_index field
  procedure :: index_i
    !! Return index_i field
  procedure :: index_j
    !! Return index_j field

  procedure, private :: set_index

#if FCKIT_FINAL_NOT_INHERITING
  final :: StructuredColumns__final_auto
#endif

END TYPE atlas_functionspace_StructuredColumns

interface atlas_functionspace_StructuredColumns
  module procedure StructuredColumns__cptr
  module procedure StructuredColumns__grid
  module procedure StructuredColumns__grid_dist
  module procedure StructuredColumns__grid_dist_levels
end interface


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

subroutine assignment_operator_hook(this,other)
  class(atlas_functionspace_StructuredColumns) :: this
  class(fckit_owned_object) :: other
  call this%set_index()
  FCKIT_SUPPRESS_UNUSED(other)
end subroutine

function StructuredColumns__cptr(cptr) result(this)
  type(atlas_functionspace_StructuredColumns) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%set_index()
  call this%return()
end function

function empty_config() result(config)
  type(atlas_Config) :: config
  config = atlas_Config()
  call config%return()
end function

function StructuredColumns__grid(grid, halo, levels) result(this)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_functionspace_StructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  integer, optional :: halo
  integer, optional :: levels
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this instead of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__functionspace__StructuredColumns__new__grid( grid%c_ptr(), config%c_ptr() ) )
  call this%set_index()
  call config%final()
  call this%return()
end function

function StructuredColumns__grid_dist(grid, distribution, halo, levels) result(this)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_functionspace_StructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_griddistribution), intent(in) :: distribution
  integer, optional :: halo
  integer, optional :: levels
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this instead of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__functionspace__StructuredColumns__new__grid_dist( &
      & grid%c_ptr(), distribution%c_ptr(), config%c_ptr() ) )
  call this%set_index()
  call config%final()
  call this%return()
end function

function StructuredColumns__grid_dist_levels(grid, distribution, levels, halo) result(this)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_functionspace_StructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_griddistribution), intent(in) :: distribution
  integer, optional :: halo
  real(c_double) :: levels(:)
  type(atlas_Config) :: config
  type(atlas_Vertical) :: vertical
  config = empty_config() ! Due to PGI compiler bug, we have to do this insted of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  call config%set("levels",size(levels))
  vertical = atlas_Vertical(levels)
  call this%reset_c_ptr( atlas__functionspace__StructuredColumns__new__grid_dist_vert( &
      & grid%c_ptr(), distribution%c_ptr(), vertical%c_ptr(), config%c_ptr() ) )
  call this%set_index()
  call config%final()
  call vertical%final()
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
  use, intrinsic :: iso_c_binding
  use atlas_functionspace_StructuredColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer(ATLAS_KIND_IDX) :: checksum_size
  integer(c_int) :: checksum_allocated
  call atlas__fs__StructuredColumns__checksum_fieldset( &
    & this%c_ptr(),fieldset%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function


function checksum_field(this,field) result(checksum)
  use, intrinsic :: iso_c_binding
  use atlas_functionspace_StructuredColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer(ATLAS_KIND_IDX) :: checksum_size
  integer(c_int) :: checksum_allocated
  call atlas__fs__StructuredColumns__checksum_field( &
      & this%c_ptr(),field%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

subroutine set_index(this)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  use atlas_functionspace_StructuredColumns_c_binding
  class(atlas_functionspace_StructuredColumns), intent(inout) :: this
  type(c_ptr) :: index_cptr
  integer(ATLAS_KIND_IDX), pointer :: index_fptr(:)
  integer(ATLAS_KIND_IDX) :: i_min, i_max, j_min, j_max
  integer(ATLAS_KIND_IDX) :: ni, nj
  call atlas__fs__StructuredColumns__index_host( this%c_ptr(), index_cptr, i_min, i_max, j_min, j_max )
  ni = i_max-i_min+1;
  nj = j_max-j_min+1;
  call c_f_pointer( index_cptr, index_fptr, (/ni*nj/) )
  this%index(i_min:i_max,j_min:j_max) => index_fptr(:)
end subroutine

function j_begin(this) result(j)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_begin(this%c_ptr())
end function

function j_end(this) result(j)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_end(this%c_ptr())
end function

function i_begin(this,j) result(i)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_begin(this%c_ptr(),j)
end function

function i_end(this,j) result(i)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_end(this%c_ptr(),j)
end function


function j_begin_halo(this) result(j)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_begin_halo(this%c_ptr())
end function

function j_end_halo(this) result(j)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  j = atlas__fs__StructuredColumns__j_end_halo(this%c_ptr())
end function

function i_begin_halo(this,j) result(i)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_begin_halo(this%c_ptr(),j)
end function

function i_end_halo(this,j) result(i)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  i = atlas__fs__StructuredColumns__i_end_halo(this%c_ptr(),j)
end function

function get_size(this) result(size)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: size
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  size = atlas__fs__StructuredColumns__size(this%c_ptr())
end function

function get_size_owned(this) result(size)
  use atlas_functionspace_StructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: size
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  size = atlas__fs__StructuredColumns__sizeOwned(this%c_ptr())
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

function index_i(this) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__StructuredColumns__index_i(this%c_ptr()) )
  call field%return()
end function

function index_j(this) result(field)
  use atlas_functionspace_StructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_StructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__StructuredColumns__index_j(this%c_ptr()) )
  call field%return()
end function

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine StructuredColumns__final_auto(this)
  type(atlas_functionspace_StructuredColumns), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_functionspace_StructuredColumns__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

end module atlas_functionspace_StructuredColumns_module

