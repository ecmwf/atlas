#include "atlas/atlas_f.h"

module atlas_functionspace_NodeColumns_module

use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Mesh_module, only: atlas_Mesh
use atlas_mesh_Nodes_module, only: atlas_mesh_Nodes
use atlas_GatherScatter_module, only: atlas_GatherScatter
use atlas_HaloExchange_module, only: atlas_HaloExchange
use atlas_Checksum_module, only: atlas_Checksum
use atlas_Config_module, only: atlas_Config
use atlas_kinds_module, only: ATLAS_KIND_GIDX

implicit none

private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_mesh_Nodes
private :: atlas_GatherScatter
private :: atlas_HaloExchange
private :: atlas_Checksum
private :: atlas_Mesh
private :: atlas_Config
private :: ATLAS_KIND_GIDX

public :: atlas_functionspace_NodeColumns

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_NodeColumns

! Purpose :
! -------
!   *atlas_functionspace_NodeColumns* : Interpretes fields defined in nodes

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains


  procedure, public :: nb_nodes
  procedure, public :: mesh
  procedure, public :: nodes

  procedure, private :: halo_exchange_fieldset
  procedure, private :: halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_fieldset, halo_exchange_field
  procedure, public :: get_halo_exchange

  procedure, private :: gather_fieldset
  procedure, private :: gather_field
  generic, public :: gather => gather_fieldset, gather_field
  procedure, public :: get_gather

  procedure, private :: scatter_fieldset
  procedure, private :: scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field
  procedure, public :: get_scatter

  procedure, private :: checksum_fieldset
  procedure, private :: checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field
  procedure, public :: get_checksum

  procedure, private :: sum_real64_r0
  procedure, private :: sum_real32_r0
  procedure, private :: sum_int64_r0
  procedure, private :: sum_int32_r0
  procedure, private :: sum_real64_r1
  procedure, private :: sum_real32_r1
  procedure, private :: sum_int64_r1
  procedure, private :: sum_int32_r1
  procedure, private :: order_independent_sum_real32_r0
  procedure, private :: order_independent_sum_real64_r0
  procedure, private :: order_independent_sum_real32_r1
  procedure, private :: order_independent_sum_real64_r1
  procedure, private :: minimum_real32_r0
  procedure, private :: minimum_real64_r0
  procedure, private :: minimum_int32_r0
  procedure, private :: minimum_int64_r0
  procedure, private :: minimum_real32_r1
  procedure, private :: minimum_real64_r1
  procedure, private :: minimum_int32_r1
  procedure, private :: minimum_int64_r1
  procedure, private :: maximum_real32_r0
  procedure, private :: maximum_real64_r0
  procedure, private :: maximum_int32_r0
  procedure, private :: maximum_int64_r0
  procedure, private :: maximum_real32_r1
  procedure, private :: maximum_real64_r1
  procedure, private :: maximum_int32_r1
  procedure, private :: maximum_int64_r1
  procedure, private :: minloc_real32_r0
  procedure, private :: minloc_real64_r0
  procedure, private :: minloc_int32_r0
  procedure, private :: minloc_int64_r0
  procedure, private :: minloc_real32_r1
  procedure, private :: minloc_real64_r1
  procedure, private :: minloc_int32_r1
  procedure, private :: minloc_int64_r1
  procedure, private :: maxloc_real32_r0
  procedure, private :: maxloc_real64_r0
  procedure, private :: maxloc_int32_r0
  procedure, private :: maxloc_int64_r0
  procedure, private :: maxloc_real32_r1
  procedure, private :: maxloc_real64_r1
  procedure, private :: maxloc_int32_r1
  procedure, private :: maxloc_int64_r1
  procedure, private :: mean_real32_r0
  procedure, private :: mean_real64_r0
  procedure, private :: mean_int32_r0
  procedure, private :: mean_int64_r0
  procedure, private :: mean_real32_r1
  procedure, private :: mean_real64_r1
  procedure, private :: mean_int32_r1
  procedure, private :: mean_int64_r1
  procedure, private :: mean_and_stddev_real32_r0
  procedure, private :: mean_and_stddev_real64_r0
  procedure, private :: mean_and_stddev_int32_r0
  procedure, private :: mean_and_stddev_int64_r0
  procedure, private :: mean_and_stddev_real32_r1
  procedure, private :: mean_and_stddev_real64_r1
  procedure, private :: mean_and_stddev_int32_r1
  procedure, private :: mean_and_stddev_int64_r1

  generic, public :: minimum => &
    & minimum_real32_r0, minimum_real32_r1, &
    & minimum_real64_r0, minimum_real64_r1, &
    & minimum_int32_r0,  minimum_int32_r1,  &
    & minimum_int64_r0,  minimum_int64_r1

  procedure, public :: minimum_per_level

  generic, public :: maximum => &
    & maximum_real32_r0, maximum_real32_r1, &
    & maximum_real64_r0, maximum_real64_r1, &
    & maximum_int32_r0,  maximum_int32_r1,  &
    & maximum_int64_r0,  maximum_int64_r1

  procedure, public :: maximum_per_level

  generic, public :: minimum_and_location => &
    & minloc_real32_r0, minloc_real32_r1, &
    & minloc_real64_r0, minloc_real64_r1, &
    & minloc_int32_r0,  minloc_int32_r1,  &
    & minloc_int64_r0,  minloc_int64_r1

  procedure, public :: minimum_and_location_per_level => &
    & minloc_per_level

  generic, public :: maximum_and_location => &
    & maxloc_real32_r0, maxloc_real32_r1, &
    & maxloc_real64_r0, maxloc_real64_r1, &
    & maxloc_int32_r0,  maxloc_int32_r1,  &
    & maxloc_int64_r0,  maxloc_int64_r1

  procedure, public :: maximum_and_location_per_level => &
    & maxloc_per_level

  generic, public :: sum => &
    & sum_real32_r0, sum_real32_r1, &
    & sum_real64_r0, sum_real64_r1, &
    & sum_int32_r0,  sum_int32_r1,  &
    & sum_int64_r0,  sum_int64_r1

  procedure, public :: sum_per_level

  generic, public :: order_independent_sum => &
    & order_independent_sum_real32_r0, order_independent_sum_real32_r1, &
    & order_independent_sum_real64_r0, order_independent_sum_real64_r1


  procedure, public :: order_independent_sum_per_level

  generic, public :: mean => &
    & mean_real32_r0, mean_real32_r1, &
    & mean_real64_r0, mean_real64_r1, &
    & mean_int32_r0,  mean_int32_r1,  &
    & mean_int64_r0,  mean_int64_r1

  procedure, public :: mean_per_level

  generic, public :: mean_and_standard_deviation => &
    & mean_and_stddev_real32_r0, mean_and_stddev_real32_r1, &
    & mean_and_stddev_real64_r0, mean_and_stddev_real64_r1, &
    & mean_and_stddev_int32_r0,  mean_and_stddev_int32_r1,  &
    & mean_and_stddev_int64_r0,  mean_and_stddev_int64_r1

  procedure, public :: mean_and_standard_deviation_per_level => &
    & mean_and_stddev_per_level

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_functionspace_NodeColumns__final_auto
#endif

END TYPE atlas_functionspace_NodeColumns

interface atlas_functionspace_NodeColumns
  module procedure constructor__cptr
  module procedure constructor
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

!------------------------------------------------------------------------------

function constructor__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_functionspace_NodeColumns) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

!------------------------------------------------------------------------------

function constructor(mesh,halo,levels) result(this)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_functionspace_NodeColumns) :: this
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  integer, intent(in), optional :: levels
  type(atlas_Config) :: config
  config = atlas_Config()
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__NodesFunctionSpace__new(mesh%c_ptr(),config%c_ptr()) )
  call config%final()
  call this%return()
end function

!------------------------------------------------------------------------------

function nb_nodes(this)
  use atlas_functionspace_NodeColumns_c_binding
  integer :: nb_nodes
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  nb_nodes = atlas__NodesFunctionSpace__nb_nodes(this%c_ptr())
end function

!------------------------------------------------------------------------------

function mesh(this)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  call mesh%reset_c_ptr( atlas__NodesFunctionSpace__mesh(this%c_ptr()) )
  call mesh%return()
end function

!------------------------------------------------------------------------------

function nodes(this)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_mesh_Nodes) :: nodes
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  call nodes%reset_c_ptr( atlas__NodesFunctionSpace__nodes(this%c_ptr()) )
  call nodes%return()
end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine halo_exchange_fieldset(this,fieldset)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__NodesFunctionSpace__halo_exchange_fieldset(this%c_ptr(),fieldset%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine halo_exchange_field(this,field)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__NodesFunctionSpace__halo_exchange_field(this%c_ptr(),field%c_ptr())
end subroutine

!------------------------------------------------------------------------------

function get_gather(this) result(gather)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_GatherScatter) :: gather
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  call gather%reset_c_ptr( atlas__NodesFunctioNSpace__get_gather(this%c_ptr()) )
end function

!------------------------------------------------------------------------------

function get_scatter(this) result(gather)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_GatherScatter) :: gather
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  call gather%reset_c_ptr( atlas__NodesFunctioNSpace__get_scatter(this%c_ptr()) )
end function

!------------------------------------------------------------------------------

subroutine gather_fieldset(this,local,global)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__NodesFunctionSpace__gather_fieldset(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine gather_field(this,local,global)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__NodesFunctionSpace__gather_field(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine scatter_fieldset(this,global,local)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__NodesFunctionSpace__scatter_fieldset(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine scatter_field(this,global,local)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__NodesFunctionSpace__scatter_field(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

!------------------------------------------------------------------------------

function get_halo_exchange(this) result(halo_exchange)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_HaloExchange) :: halo_exchange
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  call halo_exchange%reset_c_ptr( atlas__NodesFunctioNSpace__get_halo_exchange(this%c_ptr()) )
end function

!------------------------------------------------------------------------------

function get_checksum(this) result(checksum)
  use atlas_functionspace_NodeColumns_c_binding
  type(atlas_Checksum) :: checksum
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  call checksum%reset_c_ptr( atlas__NodesFunctioNSpace__get_checksum(this%c_ptr()) )
end function

!------------------------------------------------------------------------------

function checksum_fieldset(this,fieldset) result(checksum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__NodesFunctionSpace__checksum_fieldset(this%c_ptr(),fieldset%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

!------------------------------------------------------------------------------

function checksum_field(this,field) result(checksum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__NodesFunctionSpace__checksum_field(this%c_ptr(),field%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

!------------------------------------------------------------------------------

subroutine minimum_real32_r0(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: minimum
  call atlas__NodesFunctionSpace__min_float(this%c_ptr(),field%c_ptr(),minimum)
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_real32_r1(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float,c_ptr,c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: minimum(:)
  type(c_ptr) :: min_cptr
  real(c_float), pointer :: min_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__min_arr_float(this%c_ptr(),field%c_ptr(),min_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  allocate(minimum(min_size))
  minimum(:) = min_fptr(:)
  call c_ptr_free(min_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_real32_r0(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: maximum
  call atlas__NodesFunctionSpace__max_float(this%c_ptr(),field%c_ptr(),maximum)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_real32_r1(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: maximum(:)
  type(c_ptr) :: max_cptr
  real(c_float), pointer :: max_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__max_arr_float(this%c_ptr(),field%c_ptr(),max_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  allocate(maximum(max_size))
  maximum(:) = max_fptr(:)
  call c_ptr_free(max_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_real32_r0(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloc_float(this%c_ptr(),field%c_ptr(),minimum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_real32_r0(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloc_float(this%c_ptr(),field%c_ptr(),maximum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_real32_r1(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_ptr, c_f_pointer, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: min_cptr, loc_cptr
  real(c_float), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloc_arr_float(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_real32_r1(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_ptr, c_f_pointer, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: max_cptr, loc_cptr
  real(c_float), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloc_arr_float(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine sum_real32_r0(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_float(this%c_ptr(),field%c_ptr(),sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine sum_real32_r1(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_ptr_free
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_float), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_arr_float(this%c_ptr(),field%c_ptr(),sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call c_ptr_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine order_independent_sum_real32_r0(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_float(this%c_ptr(),field%c_ptr(),sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine order_independent_sum_real32_r1(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_ptr_free
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_float), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_arr_float(this%c_ptr(),field%c_ptr(),sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call c_ptr_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_real32_r0(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_float(this%c_ptr(),field%c_ptr(),mean,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_real32_r1(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_ptr_free
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: mean(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr
  real(c_float), pointer :: mean_fptr(:)
  integer :: mean_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_arr_float(this%c_ptr(),field%c_ptr(),mean_cptr,mean_size,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/mean_size/))
  allocate(mean(mean_size))
  mean(:) = mean_fptr(:)
  call c_ptr_free(mean_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_real32_r0(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: mean
  real(c_float), intent(out) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_float(this%c_ptr(),field%c_ptr(),mean,stddev,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_real32_r1(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int ,c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_ptr_free
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: mean(:)
  real(c_float), allocatable, intent(out) :: stddev(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr, stddev_cptr
  real(c_float), pointer :: mean_fptr(:), stddev_fptr(:)
  integer :: varsize
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_arr_float(this%c_ptr(),field%c_ptr(),mean_cptr,stddev_cptr,varsize,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/varsize/))
  call c_f_pointer(stddev_cptr,stddev_fptr,(/varsize/))
  allocate(mean(varsize))
  allocate(stddev(varsize))
  mean(:) = mean_fptr(:)
  stddev(:) = stddev_fptr(:)
  call c_ptr_free(mean_cptr)
  call c_ptr_free(stddev_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_real64_r0(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: minimum
  call atlas__NodesFunctionSpace__min_double(this%c_ptr(),field%c_ptr(),minimum)
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_real64_r1(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_ptr_free
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: minimum(:)
  type(c_ptr) :: min_cptr
  real(c_double), pointer :: min_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__min_arr_double(this%c_ptr(),field%c_ptr(),min_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  allocate(minimum(min_size))
  minimum(:) = min_fptr(:)
  call c_ptr_free(min_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_real64_r0(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: maximum
  call atlas__NodesFunctionSpace__max_double(this%c_ptr(),field%c_ptr(),maximum)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_real64_r1(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_ptr, c_double, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: maximum(:)
  type(c_ptr) :: max_cptr
  real(c_double), pointer :: max_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__max_arr_double(this%c_ptr(),field%c_ptr(),max_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  allocate(maximum(max_size))
  maximum(:) = max_fptr(:)
  call c_ptr_free(max_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_real64_r0(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloc_double(this%c_ptr(),field%c_ptr(),minimum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_real64_r0(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloc_double(this%c_ptr(),field%c_ptr(),maximum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_real64_r1(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_long, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: min_cptr, loc_cptr
  real(c_double), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloc_arr_double(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_real64_r1(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_ptr, c_f_pointer, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: max_cptr, loc_cptr
  real(c_double), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloc_arr_double(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine sum_real64_r0(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_double(this%c_ptr(),field%c_ptr(),sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine sum_real64_r1(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_double), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_arr_double(this%c_ptr(),field%c_ptr(),sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call c_ptr_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine order_independent_sum_real64_r0(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_double(this%c_ptr(),field%c_ptr(),sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine order_independent_sum_real64_r1(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  real(c_double), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_arr_double(this%c_ptr(),field%c_ptr(),sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call c_ptr_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_real64_r0(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_double(this%c_ptr(),field%c_ptr(),mean,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_real64_r1(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: mean(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr
  real(c_double), pointer :: mean_fptr(:)
  integer :: mean_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_arr_double(this%c_ptr(),field%c_ptr(),mean_cptr,mean_size,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/mean_size/))
  allocate(mean(mean_size))
  mean(:) = mean_fptr(:)
  call c_ptr_free(mean_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_real64_r0(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: mean
  real(c_double), intent(out) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_double(this%c_ptr(),field%c_ptr(),mean,stddev,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_real64_r1(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: mean(:)
  real(c_double), allocatable, intent(out) :: stddev(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr, stddev_cptr
  real(c_double), pointer :: mean_fptr(:), stddev_fptr(:)
  integer :: varsize
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_arr_double(this%c_ptr(),field%c_ptr(),mean_cptr,stddev_cptr,varsize,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/varsize/))
  call c_f_pointer(stddev_cptr,stddev_fptr,(/varsize/))
  allocate(mean(varsize))
  allocate(stddev(varsize))
  mean(:) = mean_fptr(:)
  stddev(:) = stddev_fptr(:)
  call c_ptr_free(mean_cptr)
  call c_ptr_free(stddev_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_int64_r0(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: minimum
  call atlas__NodesFunctionSpace__min_long(this%c_ptr(),field%c_ptr(),minimum)
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_int64_r1(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: minimum(:)
  type(c_ptr) :: min_cptr
  integer(c_long), pointer :: min_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__min_arr_long(this%c_ptr(),field%c_ptr(),min_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  allocate(minimum(min_size))
  minimum(:) = min_fptr(:)
  call c_ptr_free(min_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_int64_r0(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: maximum
  call atlas__NodesFunctionSpace__max_long(this%c_ptr(),field%c_ptr(),maximum)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_int64_r1(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: maximum(:)
  type(c_ptr) :: max_cptr
  integer(c_long), pointer :: max_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__max_arr_long(this%c_ptr(),field%c_ptr(),max_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  allocate(maximum(max_size))
  maximum(:) = max_fptr(:)
  call c_ptr_free(max_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_int64_r0(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloc_long(this%c_ptr(),field%c_ptr(),minimum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_int64_r0(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloc_long(this%c_ptr(),field%c_ptr(),maximum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_int64_r1(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer, c_long, c_ptr
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: min_cptr, loc_cptr
  integer(c_long), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloc_arr_long(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_int64_r1(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: max_cptr, loc_cptr
  integer(c_long), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloc_arr_long(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine sum_int64_r0(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_long(this%c_ptr(),field%c_ptr(),sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine sum_int64_r1(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  integer(c_long), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_arr_long(this%c_ptr(),field%c_ptr(),sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call c_ptr_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_int64_r0(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_long(this%c_ptr(),field%c_ptr(),mean,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_int64_r1(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: mean(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr
  integer(c_long), pointer :: mean_fptr(:)
  integer :: mean_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_arr_long(this%c_ptr(),field%c_ptr(),mean_cptr,mean_size,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/mean_size/))
  allocate(mean(mean_size))
  mean(:) = mean_fptr(:)
  call c_ptr_free(mean_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_int64_r0(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: mean
  integer(c_long), intent(out) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_long(this%c_ptr(),field%c_ptr(),mean,stddev,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_int64_r1(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: mean(:)
  integer(c_long), allocatable, intent(out) :: stddev(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr, stddev_cptr
  integer(c_long), pointer :: mean_fptr(:), stddev_fptr(:)
  integer :: varsize
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_arr_long(this%c_ptr(),field%c_ptr(),mean_cptr,stddev_cptr,varsize,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/varsize/))
  call c_f_pointer(stddev_cptr,stddev_fptr,(/varsize/))
  allocate(mean(varsize))
  allocate(stddev(varsize))
  mean(:) = mean_fptr(:)
  stddev(:) = stddev_fptr(:)
  call c_ptr_free(mean_cptr)
  call c_ptr_free(stddev_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_int32_r0(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: minimum
  call atlas__NodesFunctionSpace__min_int(this%c_ptr(),field%c_ptr(),minimum)
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_int32_r1(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: minimum(:)
  type(c_ptr) :: min_cptr
  integer(c_int), pointer :: min_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__min_arr_int(this%c_ptr(),field%c_ptr(),min_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  allocate(minimum(min_size))
  minimum(:) = min_fptr(:)
  call c_ptr_free(min_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_int32_r0(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: maximum
  call atlas__NodesFunctionSpace__max_int(this%c_ptr(),field%c_ptr(),maximum)
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_int32_r1(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: maximum(:)
  type(c_ptr) :: max_cptr
  integer(c_int), pointer :: max_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__max_arr_int(this%c_ptr(),field%c_ptr(),max_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  allocate(maximum(max_size))
  maximum(:) = max_fptr(:)
  call c_ptr_free(max_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_int32_r0(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloc_int(this%c_ptr(),field%c_ptr(),minimum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_int32_r0(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloc_int(this%c_ptr(),field%c_ptr(),maximum,loc)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_int32_r1(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: min_cptr, loc_cptr
  integer(c_int), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloc_arr_int(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_int32_r1(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  type(c_ptr) :: max_cptr, loc_cptr
  integer(c_int), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloc_arr_int(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine sum_int32_r0(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_int(this%c_ptr(),field%c_ptr(),sum,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine sum_int32_r1(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: sum(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: sum_cptr
  integer(c_int), pointer :: sum_fptr(:)
  integer :: sum_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_arr_int(this%c_ptr(),field%c_ptr(),sum_cptr,sum_size,opt_N)
  call c_f_pointer(sum_cptr,sum_fptr,(/sum_size/))
  allocate(sum(sum_size))
  sum(:) = sum_fptr(:)
  call c_ptr_free(sum_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_int32_r0(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_int(this%c_ptr(),field%c_ptr(),mean,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_int32_r1(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: mean(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr
  integer(c_int), pointer :: mean_fptr(:)
  integer :: mean_size
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_arr_int(this%c_ptr(),field%c_ptr(),mean_cptr,mean_size,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/mean_size/))
  allocate(mean(mean_size))
  mean(:) = mean_fptr(:)
  call c_ptr_free(mean_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_int32_r0(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: mean
  integer(c_int), intent(out) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_int(this%c_ptr(),field%c_ptr(),mean,stddev,opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_int32_r1(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: mean(:)
  integer(c_int), allocatable, intent(out) :: stddev(:)
  integer(c_int), intent(out), optional :: N
  type(c_ptr) :: mean_cptr, stddev_cptr
  integer(c_int), pointer :: mean_fptr(:), stddev_fptr(:)
  integer :: varsize
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_arr_int(this%c_ptr(),field%c_ptr(),mean_cptr,stddev_cptr,varsize,opt_N)
  call c_f_pointer(mean_cptr,mean_fptr,(/varsize/))
  call c_f_pointer(stddev_cptr,stddev_fptr,(/varsize/))
  allocate(mean(varsize))
  allocate(stddev(varsize))
  mean(:) = mean_fptr(:)
  stddev(:) = stddev_fptr(:)
  call c_ptr_free(mean_cptr)
  call c_ptr_free(stddev_cptr)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_real32_r0(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out), optional :: level
  integer(c_long) :: loc
  integer(c_int) :: opt_lev
  call atlas__NodesFunctionSpace__minloclev_float(this%c_ptr(),field%c_ptr(),minimum,loc,opt_lev)
  location = loc
  if( present(level) ) level = opt_lev
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_real32_r0(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out), optional :: level
  integer(c_long) :: loc
  integer(c_int) :: opt_lev
  call atlas__NodesFunctionSpace__maxloclev_float(this%c_ptr(),field%c_ptr(),maximum,loc,opt_lev)
  location = loc
  if( present(level) ) level = opt_lev
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_real32_r1(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_float, c_long, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out), optional :: level(:)
  type(c_ptr) :: min_cptr, loc_cptr, lev_cptr
  real(c_float), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_long),pointer :: lev_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloclev_arr_float(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,lev_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  if( present(level) ) then
    call c_f_pointer(lev_cptr,lev_fptr,(/min_size/))
    allocate(level(min_size))
    level(:) = lev_fptr(:)
  endif
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_real32_r1(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_long , c_f_pointer, c_ptr
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_float), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out), optional :: level(:)
  type(c_ptr) :: max_cptr, loc_cptr, lev_cptr
  real(c_float), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_long),pointer :: lev_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloclev_arr_float(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,lev_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  if( present(level) ) then
    call c_f_pointer(lev_cptr,lev_fptr,(/max_size/))
    allocate(level(max_size))
    level(:) = lev_fptr(:)
  endif
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_real64_r0(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_long
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out), optional :: level
  integer(c_long) :: loc
  integer(c_int) :: opt_lev
  call atlas__NodesFunctionSpace__minloclev_double(this%c_ptr(),field%c_ptr(),minimum,loc,opt_lev)
  location = loc
  if( present(level) ) level = opt_lev
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_real64_r0(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out), optional :: level
  integer(c_long) :: loc
  integer(c_int) :: opt_lev
  call atlas__NodesFunctionSpace__maxloclev_double(this%c_ptr(),field%c_ptr(),maximum,loc,opt_lev)
  location = loc
  if( present(level) ) level = opt_lev
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_real64_r1(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out), optional :: level(:)
  type(c_ptr) :: min_cptr, loc_cptr, lev_cptr
  real(c_double), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_long),pointer :: lev_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloclev_arr_double(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,lev_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  if( present(level) ) then
    call c_f_pointer(lev_cptr,lev_fptr,(/min_size/))
    allocate(level(min_size))
    level(:) = lev_fptr(:)
  endif
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_real64_r1(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_double, c_ptr, c_int, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  real(c_double), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out), optional :: level(:)
  type(c_ptr) :: max_cptr, loc_cptr, lev_cptr
  real(c_double), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_long),pointer :: lev_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloclev_arr_double(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,lev_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  if( present(level) ) then
    call c_f_pointer(lev_cptr,lev_fptr,(/max_size/))
    allocate(level(max_size))
    level(:) = lev_fptr(:)
  endif
  call c_ptr_free(loc_cptr)
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_int64_r0(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out), optional :: level
  integer(c_long) :: loc
  integer(c_int) :: opt_lev
  call atlas__NodesFunctionSpace__minloclev_long(this%c_ptr(),field%c_ptr(),minimum,loc,opt_lev)
  location = loc
  if( present(level) ) level = opt_lev
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_int64_r0(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out), optional :: level
  integer(c_long) :: loc
  integer(c_int) :: opt_lev
  call atlas__NodesFunctionSpace__maxloclev_long(this%c_ptr(),field%c_ptr(),maximum,loc,opt_lev)
  location = loc
  if( present(level) ) level = opt_lev
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_int64_r1(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int, c_ptr, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out), optional :: level(:)
  type(c_ptr) :: min_cptr, loc_cptr, lev_cptr
  integer(c_long), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_long),pointer :: lev_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloclev_arr_long(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,lev_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  if( present(level) ) then
    call c_f_pointer(lev_cptr,lev_fptr,(/min_size/))
    allocate(level(min_size))
    level(:) = lev_fptr(:)
  endif
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_int64_r1(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_long, c_int, c_ptr, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_long), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out), optional :: level(:)
  type(c_ptr) :: max_cptr, loc_cptr, lev_cptr
  integer(c_long), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_long),pointer :: lev_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloclev_arr_long(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,lev_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  if( present(level) ) then
    call c_f_pointer(lev_cptr,lev_fptr,(/max_size/))
    allocate(level(max_size))
    level(:) = lev_fptr(:)
  endif
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_int32_r0(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: minimum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out) :: level
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__minloclev_int(this%c_ptr(),field%c_ptr(),minimum,loc,level)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_int32_r0(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), intent(out) :: maximum
  integer(ATLAS_KIND_GIDX), intent(out) :: location
  integer(c_int), intent(out) :: level
  integer(c_long) :: loc
  call atlas__NodesFunctionSpace__maxloclev_int(this%c_ptr(),field%c_ptr(),maximum,loc,level)
  location = loc
end subroutine

!------------------------------------------------------------------------------

subroutine minloclev_int32_r1(this,field,minimum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: minimum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out) :: level(:)
  type(c_ptr) :: min_cptr, loc_cptr, lev_cptr
  integer(c_int), pointer :: min_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_int),pointer :: lev_fptr(:)
  integer :: min_size
  call atlas__NodesFunctionSpace__minloclev_arr_int(this%c_ptr(),field%c_ptr(),min_cptr,loc_cptr,lev_cptr,min_size)
  call c_f_pointer(min_cptr,min_fptr,(/min_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/min_size/))
  call c_f_pointer(lev_cptr,lev_fptr,(/min_size/))
  allocate(minimum(min_size))
  allocate(location(min_size))
  allocate(level(min_size))
  minimum(:) = min_fptr(:)
  location(:) = loc_fptr(:)
  level(:) = lev_fptr(:)
  call c_ptr_free(min_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine maxloclev_int32_r1(this,field,maximum,location,level)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_long, c_f_pointer
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field) :: field
  integer(c_int), allocatable, intent(out) :: maximum(:)
  integer(ATLAS_KIND_GIDX), allocatable, intent(out) :: location(:)
  integer(c_int), allocatable, intent(out) :: level(:)
  type(c_ptr) :: max_cptr, loc_cptr, lev_cptr
  integer(c_int), pointer :: max_fptr(:)
  integer(c_long),pointer :: loc_fptr(:)
  integer(c_int),pointer :: lev_fptr(:)
  integer :: max_size
  call atlas__NodesFunctionSpace__maxloclev_arr_int(this%c_ptr(),field%c_ptr(),max_cptr,loc_cptr,lev_cptr,max_size)
  call c_f_pointer(max_cptr,max_fptr,(/max_size/))
  call c_f_pointer(loc_cptr,loc_fptr,(/max_size/))
  call c_f_pointer(lev_cptr,lev_fptr,(/max_size/))
  allocate(maximum(max_size))
  allocate(location(max_size))
  allocate(level(max_size))
  maximum(:) = max_fptr(:)
  location(:) = loc_fptr(:)
  level(:) = lev_fptr(:)
  call c_ptr_free(max_cptr)
  call c_ptr_free(loc_cptr)
  call c_ptr_free(lev_cptr)
end subroutine

!------------------------------------------------------------------------------

subroutine minloc_per_level(this,field,minimum,location)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: minimum
  type(atlas_Field), intent(inout) :: location
  call atlas__NodesFunctionSpace__minloc_per_level(this%c_ptr(),field%c_ptr(),minimum%c_ptr(),location%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine maxloc_per_level(this,field,maximum,location)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: maximum
  type(atlas_Field), intent(inout) :: location
  call atlas__NodesFunctionSpace__maxloc_per_level(this%c_ptr(),field%c_ptr(),maximum%c_ptr(),location%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine minimum_per_level(this,field,minimum)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: minimum
  call atlas__NodesFunctionSpace__min_per_level(this%c_ptr(),field%c_ptr(),minimum%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine maximum_per_level(this,field,maximum)
  use atlas_functionspace_NodeColumns_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: maximum
  call atlas__NodesFunctionSpace__max_per_level(this%c_ptr(),field%c_ptr(),maximum%c_ptr())
end subroutine

!------------------------------------------------------------------------------

subroutine sum_per_level(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__sum_per_level(this%c_ptr(),field%c_ptr(),sum%c_ptr(),opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine order_independent_sum_per_level(this,field,sum,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: sum
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__oisum_per_level(this%c_ptr(),field%c_ptr(),sum%c_ptr(),opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_per_level(this,field,mean,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: mean
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_per_level(this%c_ptr(),field%c_ptr(),mean%c_ptr(),opt_N)
  if( present(N) ) N = opt_N
end subroutine

!------------------------------------------------------------------------------

subroutine mean_and_stddev_per_level(this,field,mean,stddev,N)
  use atlas_functionspace_NodeColumns_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_functionspace_NodeColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Field), intent(inout) :: mean
  type(atlas_Field), intent(inout) :: stddev
  integer(c_int), intent(out), optional :: N
  integer(c_int) :: opt_N
  call atlas__NodesFunctionSpace__mean_and_stddev_per_level( &
    & this%c_ptr(),field%c_ptr(),mean%c_ptr(),stddev%c_ptr(),opt_N)
  if( present(N) ) N = opt_N
end subroutine

!-------------------------------------------------------------------------------

subroutine atlas_functionspace_NodeColumns__final_auto(this)
  type(atlas_functionspace_NodeColumns) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_functionspace_NodeColumns__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

!------------------------------------------------------------------------------

end module atlas_functionspace_NodeColumns_module

