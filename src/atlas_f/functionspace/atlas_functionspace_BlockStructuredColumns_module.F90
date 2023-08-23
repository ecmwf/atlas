! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_functionspace_BlockStructuredColumns_module

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
use atlas_Partitioner_module, only : atlas_Partitioner
use atlas_Vertical_module, only : atlas_Vertical

implicit none

private :: c_ptr, c_double
private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Grid
private :: atlas_GridDistribution
private :: atlas_Partitioner
private :: atlas_Vertical
private :: atlas_Config
private :: fckit_owned_object

public :: atlas_functionspace_BlockStructuredColumns

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_BlockStructuredColumns

! Purpose :
! -------
!   *atlas_functionspace_BlockStructuredColumns* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015  Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
  integer(ATLAS_KIND_IDX), pointer, public :: index(:,:)

contains

  procedure, public :: assignment_operator_hook

  procedure, private :: gather_fieldset
  procedure, private :: gather_field
  generic, public :: gather => gather_fieldset, gather_field

  procedure, private :: scatter_fieldset
  procedure, private :: scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field

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
  
  procedure :: levels
  procedure :: block_begin
  procedure :: block_size
  procedure :: nproma
  procedure :: nblks

  procedure :: xy
    !! Return xy coordinate field
  procedure :: z
    !! Return z coordinate field
  procedure :: partition
    !! Return partition field
  procedure :: global_index
    !! Return global_index field
  procedure :: index_i
    !! Return index_i field
  procedure :: index_j
    !! Return index_j field

  procedure :: grid

  procedure, private :: set_index

#if FCKIT_FINAL_NOT_INHERITING
  final :: BStructuredColumns__final_auto
#endif

END TYPE atlas_functionspace_BlockStructuredColumns

interface atlas_functionspace_BlockStructuredColumns
  module procedure ctor_cptr
  module procedure ctor_grid
  module procedure ctor_grid_config
  module procedure ctor_grid_dist
  module procedure ctor_grid_dist_config
  module procedure ctor_grid_dist_levels
  module procedure ctor_grid_dist_vertical
  module procedure ctor_grid_part
  module procedure ctor_grid_part_levels
  module procedure ctor_grid_part_vertical
end interface


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

subroutine assignment_operator_hook(this,other)
  class(atlas_functionspace_BlockStructuredColumns) :: this
  class(fckit_owned_object) :: other
  call this%set_index()
  FCKIT_SUPPRESS_UNUSED(other)
end subroutine

function ctor_cptr(cptr) result(this)
  type(atlas_functionspace_BlockStructuredColumns) :: this
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

function ctor_grid(grid, halo, nproma, levels) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  integer, optional :: halo
  integer, optional :: nproma
  integer, optional :: levels
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this instead of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  if( present(nproma) )   call config%set("nproma",nproma)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid( grid%CPTR_PGIBUG_A, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call this%return()
end function

function ctor_grid_config(grid, config) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_Config), intent (in) :: config
  call this%reset_c_ptr(atlas__functionspace__BStructuredColumns__new__grid( &
      & grid%CPTR_PGIBUG_A, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call this%return()
end function

function ctor_grid_dist(grid, distribution, halo, levels) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_griddistribution), intent(in) :: distribution
  integer, optional :: halo
  integer, optional :: levels
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this instead of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid_dist( &
      & grid%CPTR_PGIBUG_A, distribution%CPTR_PGIBUG_A, config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call this%return()
end function

function ctor_grid_dist_levels(grid, distribution, levels, halo) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
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
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid_dist_vert( &
      & grid%CPTR_PGIBUG_A, distribution%CPTR_PGIBUG_A, vertical%CPTR_PGIBUG_B, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call vertical%final()
  call this%return()
end function

function ctor_grid_dist_vertical(grid, distribution, vertical, halo) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_griddistribution), intent(in) :: distribution
  integer, optional :: halo
  type(atlas_Vertical) :: vertical
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this insted of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  call config%set("levels",vertical%size())
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid_dist_vert( &
      & grid%CPTR_PGIBUG_A, distribution%CPTR_PGIBUG_A, vertical%CPTR_PGIBUG_B, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call this%return()
end function

function ctor_grid_dist_config(grid, distribution, config) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_griddistribution), intent(in) :: distribution
  type(atlas_Config), intent (in) :: config
  call this%reset_c_ptr(atlas__functionspace__BStructuredColumns__new__grid_dist_config( &
      & grid%CPTR_PGIBUG_A, distribution%CPTR_PGIBUG_A, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call this%return()
end function


function ctor_grid_part(grid, partitioner, halo, levels) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_Partitioner), intent(in) :: partitioner
  integer, optional :: halo
  integer, optional :: levels
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this instead of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid_part( &
      & grid%CPTR_PGIBUG_A, partitioner%CPTR_PGIBUG_A, config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call this%return()
end function

function ctor_grid_part_levels(grid, partitioner, levels, halo) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_Partitioner), intent(in) :: partitioner
  integer, optional :: halo
  real(c_double) :: levels(:)
  type(atlas_Config) :: config
  type(atlas_Vertical) :: vertical
  config = empty_config() ! Due to PGI compiler bug, we have to do this insted of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  call config%set("levels",size(levels))
  vertical = atlas_Vertical(levels)
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid_part_vert( &
      & grid%CPTR_PGIBUG_A, partitioner%CPTR_PGIBUG_A, vertical%CPTR_PGIBUG_B, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call vertical%final()
  call this%return()
end function

function ctor_grid_part_vertical(grid, partitioner, vertical, halo) result(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_functionspace_BlockStructuredColumns) :: this
  class(atlas_Grid), intent(in) :: grid
  type(atlas_Partitioner), intent(in) :: partitioner
  integer, optional :: halo
  type(atlas_Vertical) :: vertical
  type(atlas_Config) :: config
  config = empty_config() ! Due to PGI compiler bug, we have to do this insted of "config = atlas_Config()""
  if( present(halo) )   call config%set("halo",halo)
  call config%set("levels",vertical%size())
  call this%reset_c_ptr( atlas__functionspace__BStructuredColumns__new__grid_part_vert( &
      & grid%CPTR_PGIBUG_A, partitioner%CPTR_PGIBUG_A, vertical%CPTR_PGIBUG_B, &
      & config%CPTR_PGIBUG_B ) )
  call this%set_index()
  call config%final()
  call this%return()
end function


subroutine gather_field(this,local,global)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__functionspace__BStructuredColumns__gather_field(this%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A)
end subroutine

subroutine gather_fieldset(this,local,global)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__functionspace__BStructuredColumns__gather_fieldset(this%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A)
end subroutine

subroutine scatter_field(this,global,local)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__functionspace__BStructuredColumns__scatter_field(this%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A)
end subroutine

subroutine scatter_fieldset(this,global,local)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__functionspace__BStructuredColumns__scatter_fieldset(this%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A)
end subroutine

function checksum_fieldset(this,fieldset) result(checksum)
  use, intrinsic :: iso_c_binding
  use atlas_functionspace_BlockStructuredColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer(ATLAS_KIND_IDX) :: checksum_size
  integer(c_int) :: checksum_allocated
  call atlas__fs__BStructuredColumns__checksum_fieldset( &
    & this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function


function checksum_field(this,field) result(checksum)
  use, intrinsic :: iso_c_binding
  use atlas_functionspace_BlockStructuredColumns_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer(ATLAS_KIND_IDX) :: checksum_size
  integer(c_int) :: checksum_allocated
  call atlas__fs__BStructuredColumns__checksum_field( &
      & this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

subroutine set_index(this)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  use atlas_functionspace_BlockStructuredColumns_c_binding
  class(atlas_functionspace_BlockStructuredColumns), intent(inout) :: this
  type(c_ptr) :: index_cptr
  integer(ATLAS_KIND_IDX), pointer :: index_fptr(:)
  integer(ATLAS_KIND_IDX) :: i_min, i_max, j_min, j_max
  integer(ATLAS_KIND_IDX) :: ni, nj
  call atlas__fs__BStructuredColumns__index_host( this%CPTR_PGIBUG_A, index_cptr, i_min, i_max, j_min, j_max )
  ni = i_max-i_min+1;
  nj = j_max-j_min+1;
  call c_f_pointer( index_cptr, index_fptr, (/ni*nj/) )
  this%index(i_min:i_max,j_min:j_max) => index_fptr(:)
end subroutine

function j_begin(this) result(j)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  j = atlas__fs__BStructuredColumns__j_begin(this%CPTR_PGIBUG_A)
end function

function j_end(this) result(j)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  j = atlas__fs__BStructuredColumns__j_end(this%CPTR_PGIBUG_A)
end function

function i_begin(this,j) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__i_begin(this%CPTR_PGIBUG_A,j)
end function

function i_end(this,j) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__i_end(this%CPTR_PGIBUG_A,j)
end function


function j_begin_halo(this) result(j)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  j = atlas__fs__BStructuredColumns__j_begin_halo(this%CPTR_PGIBUG_A)
end function

function j_end_halo(this) result(j)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  j = atlas__fs__BStructuredColumns__j_end_halo(this%CPTR_PGIBUG_A)
end function

function i_begin_halo(this,j) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__i_begin_halo(this%CPTR_PGIBUG_A,j)
end function

function i_end_halo(this,j) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__i_end_halo(this%CPTR_PGIBUG_A,j)
end function

function get_size(this) result(size)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: size
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  size = atlas__fs__BStructuredColumns__size(this%CPTR_PGIBUG_A)
end function

function get_size_owned(this) result(size)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: size
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  size = atlas__fs__BStructuredColumns__sizeOwned(this%CPTR_PGIBUG_A)
end function

function levels(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: levels
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  levels = atlas__fs__BStructuredColumns__levels(this%CPTR_PGIBUG_A)
end function

function xy(this) result(field)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__BStructuredColumns__xy(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function z(this) result(field)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__BStructuredColumns__z(this%CPTR_PGIBUG_A) )
  call field%return()
end function


function partition(this) result(field)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__BStructuredColumns__partition(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function global_index(this) result(field)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__BStructuredColumns__global_index(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function index_i(this) result(field)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__BStructuredColumns__index_i(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function index_j(this) result(field)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  field = atlas_Field( atlas__fs__BStructuredColumns__index_j(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function grid(this)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  type(atlas_Grid) :: grid
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  grid = atlas_Grid( atlas__fs__BStructuredColumns__grid(this%CPTR_PGIBUG_A) )
  call grid%return()
end function

function block_begin(this,j) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__block_begin(this%CPTR_PGIBUG_A,j-1) + 1
end function

function block_size(this,j) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  integer(ATLAS_KIND_IDX), intent(in) :: j
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__block_size(this%CPTR_PGIBUG_A,j-1)
end function

function nproma(this) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__nproma(this%CPTR_PGIBUG_A)
end function

function nblks(this) result(i)
  use atlas_functionspace_BlockStructuredColumns_c_binding
  integer(ATLAS_KIND_IDX) :: i
  class(atlas_functionspace_BlockStructuredColumns), intent(in) :: this
  i = atlas__fs__BStructuredColumns__nblks(this%CPTR_PGIBUG_A)
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine BStructuredColumns__final_auto(this)
  type(atlas_functionspace_BlockStructuredColumns), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_functionspace_BlockStructuredColumns__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_functionspace_BlockStructuredColumns_module

