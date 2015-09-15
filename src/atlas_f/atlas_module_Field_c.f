! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! Field routines

function atlas_Field__cptr(cptr) result(field)
  type(atlas_Field) :: field
  type(c_ptr), intent(in) :: cptr
  field%cpp_object_ptr = cptr
  call field%attach()
  call atlas_return(field)
end function

function atlas_Field__create(params) result(field)
  type(atlas_Field) :: field
  class(atlas_Config), intent(in) :: params
  field = atlas_Field__cptr( atlas__Field__create(params%cpp_object_ptr) )
  call atlas_return(field)
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

  field = atlas_Field__cptr( atlas__Field__create(params%cpp_object_ptr) )
  call atlas_delete(params)
  call atlas_return(field)
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

  field = atlas_Field__cptr( atlas__Field__create(params%cpp_object_ptr) )
  call atlas_delete(params)
  call atlas_return(field)
end function

subroutine atlas_Field__finalize(this)
  class(atlas_Field), intent(inout) :: this
  if( c_associated(this%cpp_object_ptr) ) then
    if( this%owners() <= 0 ) then
      call atlas_abort("Cannot finalize field that has no owners")
    endif
    call this%detach()
    if( this%owners() == 0 ) then
      call atlas__Field__delete(this%cpp_object_ptr)
    endif
    this%cpp_object_ptr = c_null_ptr
  endif
end subroutine

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_Field__final(this)
  type(atlas_Field), intent(inout) :: this
  call this%finalize()
end subroutine
#endif

subroutine atlas_Field__delete(this)
  class(atlas_Field), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Field__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine

subroutine atlas_Field__reset(field_out,field_in)
  type(atlas_Field), intent(inout) :: field_out
  class(atlas_Field), intent(in) :: field_in
  if( .not. atlas_compare_equal(field_out%cpp_object_ptr,field_in%cpp_object_ptr) ) then
#ifndef FORTRAN_SUPPORTS_FINAL
    call atlas_Field__finalize(field_out)
#endif
    field_out%cpp_object_ptr = field_in%cpp_object_ptr
    if( c_associated(field_out%cpp_object_ptr) ) call field_out%attach()
  endif
end subroutine


function Field__name(this) result(field_name)
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = atlas__Field__name(this%cpp_object_ptr)
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

function Field__functionspace(this) result(functionspace)
  type(atlas_NextFunctionSpace) :: functionspace
  class(atlas_Field), intent(in) :: this
  functionspace = atlas_NextFunctionSpace(atlas__Field__functionspace(this%cpp_object_ptr))
  call atlas_return(functionspace)
end function Field__functionspace

function Field__datatype(this) result(datatype)
  class(atlas_Field), intent(in) :: this
  character(len=:), allocatable :: datatype
  type(c_ptr) :: datatype_cptr
  integer(c_int) :: datatype_size
  integer(c_int) :: datatype_allocated
  call atlas__Field__datatype(this%cpp_object_ptr,datatype_cptr,datatype_size,datatype_allocated)
  write(atlas_log%msg,*) "datatype_size = ",datatype_size; call atlas_log%error()
  allocate(character(len=datatype_size) :: datatype )
  datatype= c_to_f_string_cptr(datatype_cptr)
  if( datatype_allocated == 1 ) call atlas_free(datatype_cptr)
end function Field__datatype

function Field__size(this) result(fieldsize)
  class(atlas_Field), intent(in) :: this
  integer :: fieldsize
  fieldsize = atlas__Field__size(this%cpp_object_ptr)
end function Field__size

function Field__rank(this) result(rank)
  class(atlas_Field), intent(in) :: this
  integer :: rank
  rank = atlas__Field__rank(this%cpp_object_ptr)
end function Field__rank

function Field__bytes(this) result(bytes)
  class(atlas_Field), intent(in) :: this
  real(c_double) :: bytes
  bytes = atlas__Field__bytes(this%cpp_object_ptr)
end function Field__bytes

function Field__kind(this) result(kind)
  class(atlas_Field), intent(in) :: this
  integer :: kind
  kind = atlas__Field__kind(this%cpp_object_ptr)
end function Field__kind

function Field__levels(this) result(levels)
  class(atlas_Field), intent(in) :: this
  integer :: levels
  levels = atlas__Field__levels(this%cpp_object_ptr)
end function Field__levels

function Field__metadata(this) result(metadata)
  class(atlas_Field), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  metadata%cpp_object_ptr = atlas__Field__metadata(this%cpp_object_ptr)
end function Field__metadata


function Field__shape_array(this) result(shape)
  class(atlas_Field), intent(in) :: this
  integer, pointer :: shape(:)
  type(c_ptr) :: shape_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__shapef(this%cpp_object_ptr, shape_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_int(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_int(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_int(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_long(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_long(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_long(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_float(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_float(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_float(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_float(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
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
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data2_real64_bounds

subroutine Field__access_data3_real64_bounds(this, field, field_bounds)
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data3_real64_bounds

subroutine Field__access_data4_real64_bounds(this, field, field_bounds)
  class(atlas_Field), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_shapef_double(this%cpp_object_ptr, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data4_real64_bounds
