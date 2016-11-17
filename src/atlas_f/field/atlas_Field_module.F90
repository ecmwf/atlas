#:setvar ranks [1,2,3,4]
#:setvar dim   ['',':',':,:',':,:,:',':,:,:,:',':,:,:,:,:']
#:setvar ftypes ['integer(c_int)','integer(c_long)','real(c_float)','real(c_double)', 'logical']
#:setvar ctypes ['int','long','float','double', 'int']
#:setvar dtypes ['int32', 'int64', 'real32', 'real64', 'logical32']
#:setvar types list(zip(dtypes,ftypes,ctypes))

#:def atlas_abort(string)
atlas_abort("${string}$",atlas_code_location("atlas_Field_module.F90",${_LINE_}$))
#:enddef

module atlas_field_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_long, c_double, c_float, c_f_pointer
use atlas_c_interop, only: c_to_f_string_cptr, atlas_free, strides, view1d, c_str
use atlas_refcounted_module, only : atlas_RefCounted
use atlas_Config_module, only : atlas_Config
use atlas_Logging_module, only : atlas_log
use atlas_Error_module, only: atlas_code_location, atlas_abort, atlas_throw_outofrange
implicit none

private :: c_ptr, c_int, c_long, c_double, c_float, c_f_pointer
private :: c_to_f_string_cptr, atlas_free, c_str
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
  generic :: shape => shape_array, shape_idx

  procedure :: rename
  procedure :: set_levels
  procedure :: set_functionspace

  procedure, public :: delete => atlas_Field__delete
  procedure, public :: copy => atlas_Field__copy
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_Field__final
#endif

#:for rank in ranks
#:for dtype in dtypes
  procedure, private :: access_host_data_${dtype}$_r${rank}$
  procedure, private :: access_host_data_${dtype}$_r${rank}$_shape
  procedure, private :: access_device_data_${dtype}$_r${rank}$
  procedure, private :: access_device_data_${dtype}$_r${rank}$_shape
#:endfor
#:endfor

  generic, public :: data => &
#:for rank in ranks
#:for dtype in dtypes
      & access_host_data_${dtype}$_r${rank}$, &
      & access_host_data_${dtype}$_r${rank}$_shape, &
#:endfor
#:endfor
      & dummy

  generic, public :: host_data => &
#:for rank in ranks
#:for dtype in dtypes
      & access_host_data_${dtype}$_r${rank}$, &
      & access_host_data_${dtype}$_r${rank}$_shape, &
#:endfor
#:endfor
      & dummy

  generic, public :: device_data => &
#:for rank in ranks
#:for dtype in dtypes
      & access_device_data_${dtype}$_r${rank}$, &
      & access_device_data_${dtype}$_r${rank}$_shape, &
#:endfor
#:endfor
      & dummy

  procedure, public :: is_on_host
  procedure, public :: is_on_device
  procedure, public :: clone_to_device
  procedure, public :: clone_from_device

  procedure, private :: dummy

END TYPE atlas_Field

interface atlas_Field
  module procedure atlas_Field__cptr
  module procedure atlas_Field__create
  module procedure atlas_Field__create_name_kind_shape_int32
  module procedure atlas_Field__create_name_kind_shape_int64
  module procedure atlas_Field__create_kind_shape_int32
  module procedure atlas_Field__create_kind_shape_int64

#:for rank in ranks
#:for dtype in dtypes[:-1]  #! skip logical types
  module procedure atlas_Field__wrap_${dtype}$_r${rank}$
  module procedure atlas_Field__wrap_name_${dtype}$_r${rank}$
#:endfor
#:endfor
end interface

! ----------------------------------------------------
! ENUM DataType
integer, private, parameter :: ATLAS_KIND_INT32  = -4
integer, private, parameter :: ATLAS_KIND_INT64  = -8
integer, private, parameter :: ATLAS_KIND_REAL32 =  4
integer, private, parameter :: ATLAS_KIND_REAL64 =  8
! ----------------------------------------------------


interface array_c_to_f
#:for rank in ranks
#:for dtype in dtypes
  module procedure array_c_to_f_${dtype}$_r${rank}$
#:endfor
#:endfor
end interface
!-------------------------------------------------------------------------------


!========================================================
contains
!========================================================

#:for rank  in ranks
#:for dtype,ftype,ctype in types
!-------------------------------------------------------------------------------
subroutine array_c_to_f_${dtype}$_r${rank}$(array_cptr,rank,shape_cptr,strides_cptr,array_fptr)
  type(c_ptr), intent(in) :: array_cptr
  integer(c_int), intent(in) :: rank
  type(c_ptr), intent(in) :: shape_cptr
  type(c_ptr), intent(in) :: strides_cptr
  ${ftype}$, pointer, intent(out) :: array_fptr(${dim[rank]}$)
  ${ftype}$, pointer :: tmp(${dim[rank+1]}$)
  integer, pointer :: shape(:)
  integer, pointer :: strides(:)
  integer :: eshape(0:${rank}$)
  integer :: accumulated, factor, j

  if( rank /= ${rank}$ ) call ${atlas_abort("Rank mismatch")}$

  call c_f_pointer ( shape_cptr,   shape ,   [rank] )
  call c_f_pointer ( strides_cptr, strides , [rank] )

  eshape(0)=1
  accumulated = 1
  do j=1,rank
    accumulated = accumulated*shape(j)
    factor = shape(j)*strides(j)/max(accumulated,1)
    eshape(j-1) = eshape(j-1)*factor
    eshape(j)   = shape(j)
    accumulated = accumulated*factor
  enddo
  call c_f_pointer ( array_cptr , tmp , shape=eshape )
  #{if rank == 1}# array_fptr => tmp(1,1:shape(1)) #{endif}#
  #{if rank == 2}# array_fptr => tmp(1,1:shape(1),1:shape(2)) #{endif}#
  #{if rank == 3}# array_fptr => tmp(1,1:shape(1),1:shape(2),1:shape(3)) #{endif}#
  #{if rank == 4}# array_fptr => tmp(1,1:shape(1),1:shape(2),1:shape(3),1:shape(4)) #{endif}#
end subroutine
!-------------------------------------------------------------------------------

#:endfor
#:endfor
#:for rank  in ranks
#:for dtype,ftype,ctype in types
subroutine access_host_data_${dtype}$_r${rank}$(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  ${ftype}$, pointer, intent(out) :: field(${dim[rank]}$)
  type(c_ptr) :: field_cptr
  type(c_ptr) :: shape_cptr
  type(c_ptr) :: strides_cptr
  integer(c_int) :: rank
  call atlas__Field__host_data_${ctype}$_specf(this%c_ptr(), field_cptr, rank, shape_cptr, strides_cptr)
  call array_c_to_f(field_cptr,rank,shape_cptr,strides_cptr, field)
end subroutine

subroutine access_device_data_${dtype}$_r${rank}$(this, field)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  ${ftype}$, pointer, intent(out) :: field(${dim[rank]}$)
  type(c_ptr) :: field_cptr
  type(c_ptr) :: shape_cptr
  type(c_ptr) :: strides_cptr
  integer(c_int) :: rank
  call atlas__Field__device_data_${ctype}$_specf(this%c_ptr(), field_cptr, rank, shape_cptr, strides_cptr)
  call array_c_to_f(field_cptr,rank,shape_cptr,strides_cptr, field)
end subroutine

!-------------------------------------------------------------------------------

#:endfor
#:endfor
#:for rank  in ranks
#:for dtype,ftype,ctype in types
subroutine access_host_data_${dtype}$_r${rank}$_shape(this, field, shape)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  ${ftype}$, pointer, intent(out) :: field(${dim[rank]}$)
  integer(c_int), intent(in) :: shape(:)
  type(c_ptr) :: field_cptr
  type(c_ptr) :: shape_cptr
  type(c_ptr) :: strides_cptr
  integer(c_int) :: rank
  call atlas__Field__host_data_${ctype}$_specf(this%c_ptr(), field_cptr, rank, shape_cptr, strides_cptr)
  call c_f_pointer( field_cptr, field, shape )
end subroutine

subroutine access_device_data_${dtype}$_r${rank}$_shape(this, field, shape)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  ${ftype}$, pointer, intent(out) :: field(${dim[rank]}$)
  integer(c_int), intent(in) :: shape(:)
  type(c_ptr) :: field_cptr
  type(c_ptr) :: shape_cptr
  type(c_ptr) :: strides_cptr
  integer(c_int) :: rank
  call atlas__Field__device_data_${ctype}$_specf(this%c_ptr(), field_cptr, rank, shape_cptr, strides_cptr)
  call c_f_pointer( field_cptr, field, shape )
end subroutine

!-------------------------------------------------------------------------------
#:endfor
#:endfor
subroutine dummy(this)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
end subroutine

!-------------------------------------------------------------------------------

integer function atlas_real(kind)
  integer :: kind
  if (kind == c_double) then
    atlas_real = ATLAS_KIND_REAL64
  else if (kind == c_float) then
    atlas_real = ATLAS_KIND_REAL32
  else
    call ${atlas_abort("Unsupported real kind")}$
  end if
end function

!-------------------------------------------------------------------------------

integer function atlas_integer(kind)
  integer, optional :: kind
  atlas_integer = ATLAS_KIND_INT32
  if ( present(kind) ) then
    if (kind == c_int) then
      atlas_integer = ATLAS_KIND_INT32
    else if (kind == c_long) then
      atlas_integer = ATLAS_KIND_INT64
    else
      call ${atlas_abort("Unsupported real kind")}$
    end if
  end if
end function

!-------------------------------------------------------------------------------

integer function atlas_logical(kind)
  integer, optional :: kind
  atlas_logical = ATLAS_KIND_INT32
end function

!-------------------------------------------------------------------------------

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
    call ${atlas_abort("cannot convert kind to data_type")}$
  endif
end function

!-------------------------------------------------------------------------------

function atlas_Field__cptr(cptr) result(field)
  type(atlas_Field) :: field
  type(c_ptr), intent(in) :: cptr
  call field%reset_c_ptr( cptr )
end function

!-------------------------------------------------------------------------------

function atlas_Field__create(params) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  class(atlas_Config), intent(in) :: params
  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call field%return()
end function

!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------

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


!-------------------------------------------------------------------------------

#:for rank in ranks
#:for dtype, ftype, ctype in types[:-1]  #! skip logical types
function atlas_Field__wrap_name_${dtype}$_r${rank}$(name,data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  ${ftype}$, intent(in) :: data(${dim[rank]}$)
  integer(c_int) :: shapef(${rank}$)
  integer(c_int) :: stridesf(${rank}$)
  ${ftype}$, pointer :: data1d(:)
  shapef = shape(data)
  stridesf = strides(data)
  data1d => view1d(data)
  field = atlas_Field__cptr( atlas__Field__wrap_${ctype}$_specf(name,data1d,size(shapef),shapef, stridesf) )
  call field%return()
end function
function atlas_Field__wrap_${dtype}$_r${rank}$(data) result(field)
  use atlas_field_c_binding
  type(atlas_Field) :: field
  ${ftype}$, intent(in) :: data(${dim[rank]}$)
  integer(c_int) :: shapef(${rank}$)
  integer(c_int) :: stridesf(${rank}$)
  ${ftype}$, pointer :: data1d(:)
  shapef = shape(data)
  stridesf = strides(data)
  data1d => view1d(data)
  field = atlas_Field__cptr( atlas__Field__wrap_${ctype}$_specf(c_str(""),data1d,size(shapef),shapef, stridesf) )
  call field%return()
end function
#:endfor
#:endfor

!-------------------------------------------------------------------------------

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_Field__final(this)
  type(atlas_Field), intent(inout) :: this
  call this%final()
end subroutine
#endif

!-------------------------------------------------------------------------------

subroutine atlas_Field__delete(this)
  use atlas_field_c_binding
  class(atlas_Field), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Field__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

!-------------------------------------------------------------------------------

subroutine atlas_Field__copy(this,obj_in)
  class(atlas_Field), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine

!-------------------------------------------------------------------------------

function Field__name(this) result(field_name)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = atlas__Field__name(this%c_ptr())
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

!-------------------------------------------------------------------------------

function Field__functionspace(this) result(functionspace)
  use atlas_field_c_binding
  use atlas_functionspace_module
  type(atlas_FunctionSpace) :: functionspace
  class(atlas_Field), intent(in) :: this
  functionspace = atlas_FunctionSpace(atlas__Field__functionspace(this%c_ptr()))
  call functionspace%return()
end function Field__functionspace

!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------

function Field__size(this) result(fieldsize)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: fieldsize
  fieldsize = atlas__Field__size(this%c_ptr())
end function Field__size

!-------------------------------------------------------------------------------

function Field__rank(this) result(rank)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: rank
  rank = atlas__Field__rank(this%c_ptr())
end function Field__rank

!-------------------------------------------------------------------------------

function Field__bytes(this) result(bytes)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  real(c_double) :: bytes
  bytes = atlas__Field__bytes(this%c_ptr())
end function Field__bytes

!-------------------------------------------------------------------------------

function Field__kind(this) result(kind)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: kind
  kind = atlas__Field__kind(this%c_ptr())
end function Field__kind

!-------------------------------------------------------------------------------

function Field__levels(this) result(levels)
  use atlas_field_c_binding
  class(atlas_Field), intent(in) :: this
  integer :: levels
  levels = atlas__Field__levels(this%c_ptr())
end function Field__levels

!-------------------------------------------------------------------------------

function Field__metadata(this) result(metadata)
  use atlas_field_c_binding
  use atlas_metadata_module
  class(atlas_Field), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  call metadata%reset_c_ptr( atlas__Field__metadata(this%c_ptr()) )
end function Field__metadata

!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------

subroutine set_levels(this,nb_levels)
  use atlas_field_c_binding
  class(atlas_Field), intent(inout) :: this
  integer, intent(in) :: nb_levels
  call atlas__field__set_levels(this%c_ptr(),nb_levels)
end subroutine

!-------------------------------------------------------------------------------

subroutine rename(this,name)
  use atlas_field_c_binding
  class(atlas_Field), intent(inout) :: this
  character(len=*), intent(in) :: name
  call atlas__field__rename(this%c_ptr(),c_str(name))
end subroutine

!-------------------------------------------------------------------------------

subroutine set_functionspace(this,functionspace)
  use atlas_field_c_binding
  use atlas_functionspace_module
  class(atlas_Field), intent(inout) :: this
  class(atlas_FunctionSpace), intent(in) :: functionspace
  call atlas__field__set_functionspace(this%c_ptr(),functionspace%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

function is_on_host(this)
  use atlas_field_c_binding
  logical :: is_on_host
  class(atlas_Field), intent(in) :: this
  if( atlas__Field__is_on_host(this%c_ptr()) == 1 ) then
    is_on_host = .true.
  else
    is_on_host = .false.
  endif
end function

!-------------------------------------------------------------------------------

function is_on_device(this)
  use atlas_field_c_binding
  logical :: is_on_device
  class(atlas_Field), intent(in) :: this
  if( atlas__Field__is_on_device(this%c_ptr()) == 1 ) then
    is_on_device = .true.
  else
    is_on_device = .false.
  endif
end function

!-------------------------------------------------------------------------------

subroutine clone_to_device(this)
  use atlas_field_c_binding
  class(atlas_Field), intent(inout) :: this
  call atlas__Field__clone_to_device(this%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

subroutine clone_from_device(this)
  use atlas_field_c_binding
  class(atlas_Field), intent(inout) :: this
  call atlas__Field__clone_from_device(this%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

end module atlas_field_module

