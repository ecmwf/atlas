! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! Field routines

function atlas_Field__cptr(cptr) result(field)
  type(atlas_Field) :: field
  type(c_ptr), intent(in) :: cptr
  call field%reset_c_ptr( cptr )
end function

function atlas_Field__create(params) result(field)
  type(atlas_Field) :: field
  class(atlas_Config), intent(in) :: params
  field = atlas_Field__cptr( atlas__Field__create(params%c_ptr()) )
  call field%return()
end function

function atlas_Field__create_name_kind_shape(name,kind,shape) result(field)
  type(atlas_Field) :: field
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  integer, intent(in) :: shape(:)

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

function atlas_Field__create_kind_shape(kind,shape) result(field)
  type(atlas_Field) :: field
  integer, intent(in) :: kind
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

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_Field__final(this)
  type(atlas_Field), intent(inout) :: this
  call this%final()
end subroutine
#endif

subroutine atlas_Field__delete(this)
  class(atlas_Field), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Field__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine


function Field__name(this) result(field_name)
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = atlas__Field__name(this%c_ptr())
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

function Field__functionspace(this) result(functionspace)
  type(atlas_NextFunctionSpace) :: functionspace
  class(atlas_Field), intent(in) :: this
  functionspace = atlas_NextFunctionSpace(atlas__Field__functionspace(this%c_ptr()))
  call functionspace%return()
end function Field__functionspace

function Field__datatype(this) result(datatype)
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: datatype
  type(c_ptr) :: datatype_cptr
  integer(c_int) :: datatype_size
  integer(c_int) :: datatype_allocated
  call atlas__Field__datatype(this%c_ptr(),datatype_cptr,datatype_size,datatype_allocated)
  write(atlas_log%msg,*) "datatype_size = ",datatype_size; call atlas_log%error()
  allocate(character(len=datatype_size) :: datatype )
  datatype= c_to_f_string_cptr(datatype_cptr)
  if( datatype_allocated == 1 ) call atlas_free(datatype_cptr)
end function Field__datatype

function Field__size(this) result(fieldsize)
  class(atlas_Field), intent(in) :: this
  integer :: fieldsize
  fieldsize = atlas__Field__size(this%c_ptr())
end function Field__size

function Field__rank(this) result(rank)
  class(atlas_Field), intent(in) :: this
  integer :: rank
  rank = atlas__Field__rank(this%c_ptr())
end function Field__rank

function Field__bytes(this) result(bytes)
  class(atlas_Field), intent(in) :: this
  real(c_double) :: bytes
  bytes = atlas__Field__bytes(this%c_ptr())
end function Field__bytes

function Field__kind(this) result(kind)
  class(atlas_Field), intent(in) :: this
  integer :: kind
  kind = atlas__Field__kind(this%c_ptr())
end function Field__kind

function Field__levels(this) result(levels)
  class(atlas_Field), intent(in) :: this
  integer :: levels
  levels = atlas__Field__levels(this%c_ptr())
end function Field__levels

function Field__metadata(this) result(metadata)
  class(atlas_Field), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  call metadata%reset_c_ptr( atlas__Field__metadata(this%c_ptr()) )
end function Field__metadata


function Field__shape_array(this) result(shape)
  class(atlas_Field), intent(in) :: this
  integer, pointer :: shape(:)
  type(c_ptr) :: shape_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__shapef(this%c_ptr(), shape_c_ptr, field_rank)
  call C_F_POINTER ( shape_c_ptr , shape , (/field_rank/) )
end function Field__shape_array

function Field__shape_idx(this,idx) result(shape_val)
  integer :: shape_val
  class(atlas_Field), intent(in) :: this
  integer, intent(in) :: idx
  integer, pointer :: shape(:)
  shape => this%shape_array()
  if( idx > size(shape) ) call atlas_throw_outofrange("shape",idx,size(shape),atlas_code_location(__FILE__,__LINE__))
  shape_val = shape(idx)
end function Field__shape_idx

subroutine Field__access_data1_int32(this, field)
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
  class(atlas_Field), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_float(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  if( field_rank /= 4 ) call atlas_abort("data is not of rank 4",atlas_code_location(__FILE__,__LINE__))
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real32

subroutine Field__access_data1_real64(this, field)
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
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  if( field_rank /= 4 ) call atlas_abort("data is not of rank 4",atlas_code_location(__FILE__,__LINE__))
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real64

subroutine Field__access_data2_real64_bounds(this, field, field_bounds)
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
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%c_ptr(), field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real64_bounds
