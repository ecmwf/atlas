
module atlas_field_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_long, c_double, c_float, c_f_pointer
use atlas_c_interop, only: c_to_f_string_cptr, atlas_free
use atlas_refcounted_module, only : atlas_RefCounted
use atlas_Config_module, only : atlas_Config
use atlas_Logging_module, only : atlas_log
use atlas_Error_module, only: atlas_code_location, atlas_abort, atlas_throw_outofrange
implicit none

private :: c_ptr, c_int, c_long, c_double, c_float, c_f_pointer
private :: c_to_f_string_cptr, atlas_free
private :: atlas_RefCounted, atlas_Config, atlas_log, atlas_code_location, atlas_abort, atlas_throw_outofrange

public :: atlas_Field
public :: atlas_real
public :: atlas_integer
public :: atlas_logical
public :: atlas_data_type

character(len=*), parameter :: filename = 'atlas_Field_module.F90'

private


!------------------------------------------------------------------------------
TYPE, extends(atlas_refcounted) :: atlas_Field

! Purpose :
! -------
!   *Field* : Object containing field data and Metadata

! Methods :
! -------
!   name : The name or tag this field was created with
!   data : Return the field as a fortran array of specified shape
!   Metadata : Return object that can contain a variety of Metadata

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: name => Field__name
  procedure :: functionspace => Field__functionspace
  procedure :: datatype => Field__datatype
  procedure :: metadata => Field__metadata
  procedure, private :: shape_array => Field__shape_array
  procedure, private :: shape_idx   => Field__shape_idx
  procedure :: size => Field__size
  procedure :: rank => Field__rank
  procedure :: bytes => Field__bytes
  procedure :: levels => Field__levels
  procedure :: kind => Field__kind
  procedure, private :: access_data1_logical32 => Field__access_data1_logical32
  procedure, private :: access_data2_logical32 => Field__access_data2_logical32
  procedure, private :: access_data3_logical32 => Field__access_data3_logical32
  procedure, private :: access_data1_int32 => Field__access_data1_int32
  procedure, private :: access_data2_int32 => Field__access_data2_int32
  procedure, private :: access_data3_int32 => Field__access_data3_int32
  procedure, private :: access_data1_int64 => Field__access_data1_int64
  procedure, private :: access_data2_int64 => Field__access_data2_int64
  procedure, private :: access_data3_int64 => Field__access_data3_int64
  procedure, private :: access_data1_real32 => Field__access_data1_real32
  procedure, private :: access_data2_real32 => Field__access_data2_real32
  procedure, private :: access_data3_real32 => Field__access_data3_real32
  procedure, private :: access_data4_real32 => Field__access_data4_real32
  procedure, private :: access_data1_real64 => Field__access_data1_real64
  procedure, private :: access_data2_real64 => Field__access_data2_real64
  procedure, private :: access_data3_real64 => Field__access_data3_real64
  procedure, private :: access_data4_real64 => Field__access_data4_real64
  procedure, private :: access_data2_real64_bounds => Field__access_data2_real64_bounds
  procedure, private :: access_data3_real64_bounds => Field__access_data3_real64_bounds
  procedure, private :: access_data4_real64_bounds => Field__access_data4_real64_bounds
  generic :: shape => shape_array, shape_idx
  generic :: data => &
    & access_data1_logical32, &
    & access_data1_int32, &
    & access_data1_int64, &
    & access_data1_real32, &
    & access_data1_real64, &
    & access_data2_logical32, &
    & access_data2_int32, &
    & access_data2_int64, &
    & access_data2_real32, &
    & access_data2_real64, &
    & access_data2_real64_bounds, &
    & access_data3_logical32, &
    & access_data3_int32, &
    & access_data3_int64, &
    & access_data3_real32, &
    & access_data3_real64, &
    & access_data3_real64_bounds, &
    & access_data4_real32, &
    & access_data4_real64, &
    & access_data4_real64_bounds
  procedure, public :: delete => atlas_Field__delete
  procedure, public :: copy => atlas_Field__copy
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_Field__final
#endif
END TYPE atlas_Field

interface atlas_Field
  module procedure atlas_Field__cptr
  module procedure atlas_Field__create
  module procedure atlas_Field__create_name_kind_shape_int32
  module procedure atlas_Field__create_name_kind_shape_int64
  module procedure atlas_Field__create_kind_shape_int32
  module procedure atlas_Field__create_kind_shape_int64
  module procedure atlas_Field__wrap_name_kind_shape_int32_r1
  module procedure atlas_Field__wrap_name_kind_shape_int32_r2
  module procedure atlas_Field__wrap_name_kind_shape_int32_r3
  module procedure atlas_Field__wrap_name_kind_shape_int32_r4
  module procedure atlas_Field__wrap_name_kind_shape_int64_r1
  module procedure atlas_Field__wrap_name_kind_shape_int64_r2
  module procedure atlas_Field__wrap_name_kind_shape_int64_r3
  module procedure atlas_Field__wrap_name_kind_shape_int64_r4
  module procedure atlas_Field__wrap_name_kind_shape_real32_r1
  module procedure atlas_Field__wrap_name_kind_shape_real32_r2
  module procedure atlas_Field__wrap_name_kind_shape_real32_r3
  module procedure atlas_Field__wrap_name_kind_shape_real32_r4
  module procedure atlas_Field__wrap_name_kind_shape_real64_r1
  module procedure atlas_Field__wrap_name_kind_shape_real64_r2
  module procedure atlas_Field__wrap_name_kind_shape_real64_r3
  module procedure atlas_Field__wrap_name_kind_shape_real64_r4
end interface

! ----------------------------------------------------
! ENUM DataType
integer, private, parameter :: ATLAS_KIND_INT32  = -4
integer, private, parameter :: ATLAS_KIND_INT64  = -8
integer, private, parameter :: ATLAS_KIND_REAL32 =  4
integer, private, parameter :: ATLAS_KIND_REAL64 =  8
! ----------------------------------------------------


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================



integer function atlas_real(kind)
  integer :: kind
  if (kind == c_double) then
    atlas_real = ATLAS_KIND_REAL64
  else if (kind == c_float) then
    atlas_real = ATLAS_KIND_REAL32
  else
    call atlas_abort("Unsupported real kind")
  end if
end function

integer function atlas_integer(kind)
  integer, optional :: kind
  atlas_integer = ATLAS_KIND_INT32
  if ( present(kind) ) then
    if (kind == c_int) then
      atlas_integer = ATLAS_KIND_INT32
    else if (kind == c_long) then
      atlas_integer = ATLAS_KIND_INT64
    else
      call atlas_abort("Unsupported real kind")
    end if
  end if
end function

integer function atlas_logical(kind)
  integer, optional :: kind
  atlas_logical = ATLAS_KIND_INT32
end function

function atlas_data_type(kind)
  character(len=6) :: atlas_data_type
  integer, intent(in) :: kind
  if( kind == ATLAS_KIND_INT32 ) then
    atlas_data_type = "int32"
  else if( kind == ATLAS_KIND_INT64 ) then
    atlas_data_type = "int64"
  else if( kind == ATLAS_KIND_REAL32 ) then
    atlas_data_type = "real32"
  else if( kind == ATLAS_KIND_REAL64 ) then
    atlas_data_type = "real64"
  else
    call atlas_abort("cannot convert kind to data_type", &
& atlas_code_location(filename,__LINE__))
  endif
end function



! -----------------------------------------------------------------------------
! Field routines

function atlas_Field__cptr(cptr) result(field)
  type(atlas_Field) :: field
  type(c_ptr), intent(in) :: cptr
  call field%reset_c_ptr( cptr )
end function

function atlas_Field__create(params) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  class(atlas_Config), intent(in) :: params
  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call field%return()
end function

function atlas_Field__create_name_kind_shape_int32(name,kind,shape) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  integer(c_int), intent(in) :: shape(:)

  type(atlas_Config) :: params

  params = atlas_Config()
  call params%set("creator","ArraySpec")
  call params%set("shape",shape)
  call params%set("fortran",.True.)
  call params%set("datatype",atlas_data_type(kind))
  call params%set("name",name)

  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call params%final()
  call field%return()
end function

function atlas_Field__create_name_kind_shape_int64(name,kind,shape) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  integer(c_long), intent(in) :: shape(:)

  type(atlas_Config) :: params

  params = atlas_Config()
  call params%set("creator","ArraySpec")
  call params%set("shape",shape)
  call params%set("fortran",.True.)
  call params%set("datatype",atlas_data_type(kind))
  call params%set("name",name)

  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call params%final()
  call field%return()
end function

function atlas_Field__create_kind_shape_int32(kind,shape) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  integer(c_int), intent(in) :: kind
  integer, intent(in) :: shape(:)

  type(atlas_Config) :: params

  params = atlas_Config()
  call params%set("creator","ArraySpec")
  call params%set("shape",shape)
  call params%set("fortran",.True.)
  call params%set("datatype",atlas_data_type(kind))

  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call params%final()
  call field%return()
end function

function atlas_Field__create_kind_shape_int64(kind,shape) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  integer, intent(in) :: kind
  integer(c_long), intent(in) :: shape(:)

  type(atlas_Config) :: params

  params = atlas_Config()
  call params%set("creator","ArraySpec")
  call params%set("shape",shape)
  call params%set("fortran",.True.)
  call params%set("datatype",atlas_data_type(kind))

  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call params%final()
  call field%return()
end function


function atlas_Field__wrap_name_kind_shape_int32_r1(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: data(:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(1))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_int_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int32_r2(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: data(:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(2))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_int_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int32_r3(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: data(:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(3))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_int_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int32_r4(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: data(:,:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(4))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_int_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int64_r1(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: data(:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(1))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_long_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int64_r2(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: data(:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(2))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_long_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int64_r3(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: data(:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(3))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_long_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_int64_r4(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: data(:,:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(4))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_long_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function


function atlas_Field__wrap_name_kind_shape_real32_r1(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: data(:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(1))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_float_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real32_r2(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: data(:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(2))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_float_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real32_r3(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: data(:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(3))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_float_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real32_r4(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: data(:,:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(4))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_float_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real64_r1(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: data(:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(1))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_double_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real64_r2(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: data(:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(2))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_double_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real64_r3(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: data(:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(3))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_double_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function

function atlas_Field__wrap_name_kind_shape_real64_r4(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: data(:,:,:,:)
  integer(c_int), allocatable :: shapef(:)
  allocate(shapef(4))
  shapef = shape(data)
  field = atlas_Field__cptr( atlas__Field__wrap_double_shapef(name,data,size(shapef),shapef) )
  call field%return()
end function



#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_Field__final(this)
  type(atlas_Field), intent(inout) :: this
  call this%final()
end subroutine
#endif

subroutine atlas_Field__delete(this)
  use atlas_field_c_binding
  class(atlas_Field), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Field__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_Field__copy(this,obj_in)
  class(atlas_Field), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine

function Field__name(this) result(field_name)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = atlas__Field__name(this%c_ptr())
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

function Field__functionspace(this) result(functionspace)
  use atlas_field_c_binding
  use atlas_functionspace_module
  type(atlas_FunctionSpace) :: functionspace
  class(atlas_Field), intent(in) :: this
  functionspace = atlas_FunctionSpace(atlas__Field__functionspace(this%c_ptr()))
  call functionspace%return()
end function Field__functionspace

function Field__datatype(this) result(datatype)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: datatype
  type(c_ptr) :: datatype_cptr
  integer(c_int) :: datatype_size
  integer(c_int) :: datatype_allocated
  call atlas__Field__datatype(this%c_ptr(),datatype_cptr,datatype_size,datatype_allocated)
  allocate(character(len=datatype_size) :: datatype )
  datatype= c_to_f_string_cptr(datatype_cptr)
  if( datatype_allocated == 1 ) call atlas_free(datatype_cptr)
end function Field__datatype

function Field__size(this) result(fieldsize)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: fieldsize
  fieldsize = atlas__Field__size(this%c_ptr())
end function Field__size

function Field__rank(this) result(rank)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: rank
  rank = atlas__Field__rank(this%c_ptr())
end function Field__rank

function Field__bytes(this) result(bytes)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double) :: bytes
  bytes = atlas__Field__bytes(this%c_ptr())
end function Field__bytes

function Field__kind(this) result(kind)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: kind
  kind = atlas__Field__kind(this%c_ptr())
end function Field__kind

function Field__levels(this) result(levels)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: levels
  levels = atlas__Field__levels(this%c_ptr())
end function Field__levels

function Field__metadata(this) result(metadata)
  use atlas_field_c_binding
  use atlas_metadata_module
  class(atlas_Field), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  call metadata%reset_c_ptr( atlas__Field__metadata(this%c_ptr()) )
end function Field__metadata


function Field__shape_array(this) result(shape)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer, allocatable :: shape(:)
  type(c_ptr) :: shape_c_ptr
  integer, pointer :: shape_f_ptr(:)
  integer(c_int) :: field_rank
  call atlas__Field__shapef(this%c_ptr(), shape_c_ptr, field_rank)
  call C_F_POINTER ( shape_c_ptr , shape_f_ptr , (/field_rank/) )
  allocate( shape(field_rank) )
  shape(:) = shape_f_ptr(:)
end function Field__shape_array

function Field__shape_idx(this,idx) result(shape_val)
  use atlas_field_c_binding
  integer :: shape_val
  class(atlas_Field), intent(in) :: this
  integer, intent(in) :: idx
  type(c_ptr) :: shape_c_ptr
  integer, pointer :: shape_f_ptr(:)
  integer(c_int) :: field_rank
  call atlas__Field__shapef(this%c_ptr(), shape_c_ptr, field_rank)
  call C_F_POINTER ( shape_c_ptr , shape_f_ptr , (/field_rank/) )
  if( idx > field_rank ) call atlas_throw_outofrange("shape",idx,field_rank, &
& atlas_code_location(filename,__LINE__))
  shape_val = shape_f_ptr(idx)
end function Field__shape_idx

subroutine Field__access_data1_logical32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  logical, pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_shapef_int(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_logical32

subroutine Field__access_data2_logical32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  logical, pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_int(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_logical32

subroutine Field__access_data3_logical32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  logical, pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_int(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_logical32

!--
subroutine Field__access_data1_int32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer, pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_shapef_int(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_int32

subroutine Field__access_data2_int32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer, pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_int(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_int32

subroutine Field__access_data3_int32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer, pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_int(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_int32

subroutine Field__access_data1_int64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer(c_long), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_shapef_long(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_int64

subroutine Field__access_data2_int64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer(c_long), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_long(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_int64

subroutine Field__access_data3_int64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer(c_long), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_long(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_int64

subroutine Field__access_data1_real32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_shapef_float(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_real32

subroutine Field__access_data2_real32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_float(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_real32

subroutine Field__access_data3_real32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_float(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_real32

subroutine Field__access_data4_real32(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_float(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  if( field_rank /= 4 ) call atlas_abort("data is not of rank 4", &
& atlas_code_location(filename,__LINE__))
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real32

subroutine Field__access_data1_real64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_real64

subroutine Field__access_data2_real64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_real64

subroutine Field__access_data3_real64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_real64

subroutine Field__access_data4_real64(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  if( field_rank /= 4 ) call atlas_abort("data is not of rank 4", &
& atlas_code_location(filename,__LINE__))
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real64

subroutine Field__access_data2_real64_bounds(this, field, field_bounds)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data2_real64_bounds

subroutine Field__access_data3_real64_bounds(this, field, field_bounds)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data3_real64_bounds

subroutine Field__access_data4_real64_bounds(this, field, field_bounds)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real64_bounds

end module atlas_field_module

