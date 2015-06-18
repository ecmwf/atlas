! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! Metadata routines

function atlas_Metadata__ctor() result(metadata)
  type(atlas_Metadata) :: metadata
  metadata%cpp_object_ptr = atlas__Metadata__new()
end function atlas_Metadata__ctor

subroutine atlas_Metadata__delete(this)
  class(atlas_Metadata), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Metadata__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = C_NULL_ptr
end subroutine atlas_Metadata__delete

function Metadata__has(this, name) result(value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical :: value
  integer :: value_int
  value_int =  atlas__Metadata__has(this%cpp_object_ptr, c_str(name) )
  if( value_int == 1 ) then
    value = .True.
  else
    value = .False.
  end if
end function Metadata__has

subroutine Metadata__set_logical(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical, intent(in) :: value
  integer :: value_int
  if( value ) then
    value_int = 1
  else
    value_int = 0
  end if
  call atlas__Metadata__set_int(this%cpp_object_ptr, c_str(name), value_int )
end subroutine Metadata__set_logical

subroutine Metadata__set_int32(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call atlas__Metadata__set_int(this%cpp_object_ptr, c_str(name), value)
end subroutine Metadata__set_int32

subroutine Metadata__set_real32(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Metadata__set_float(this%cpp_object_ptr, c_str(name) ,value)
end subroutine Metadata__set_real32

subroutine Metadata__set_real64(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Metadata__set_double(this%cpp_object_ptr, c_str(name) ,value)
end subroutine Metadata__set_real64

subroutine Metadata__set_string(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Metadata__set_string(this%cpp_object_ptr, c_str(name) , c_str(value) )
end subroutine Metadata__set_string

subroutine Metadata__set_mesh(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(atlas_Mesh), intent(in) :: value
  call atlas__Metadata__set_mesh(this%cpp_object_ptr, c_str(name), value%cpp_object_ptr)
end subroutine Metadata__set_mesh

subroutine Metadata__set_grid(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(atlas_ReducedGrid), intent(in) :: value
  call atlas__Metadata__set_grid(this%cpp_object_ptr, c_str(name), value%cpp_object_ptr)
end subroutine Metadata__set_grid

subroutine Metadata__set_function_space(this, name, value)
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(atlas_FunctionSpace), intent(in) :: value
  call atlas__Metadata__set_function_space(this%cpp_object_ptr, c_str(name), value%cpp_object_ptr)
end subroutine Metadata__set_function_space

subroutine Metadata__get_logical(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(out) :: value
  integer :: value_int
  value_int = atlas__Metadata__get_int(this%cpp_object_ptr,c_str(name) )
  if (value_int > 0) then
    value = .True.
  else
    value = .False.
  end if
end subroutine Metadata__get_logical

subroutine Metadata__get_int32(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(out) :: value
  value = atlas__Metadata__get_int(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_int32

subroutine Metadata__get_real32(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(out) :: value
  value = atlas__Metadata__get_float(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_real32

subroutine Metadata__get_real64(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(out) :: value
  value = atlas__Metadata__get_double(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_real64

subroutine Metadata__get_string(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(out) :: value
  character(len=MAX_STR_LEN) :: value_cstr
  call atlas__Metadata__get_string(this%cpp_object_ptr, c_str(name), value_cstr, MAX_STR_LEN )
  value = c_to_f_string_str(value_cstr)
end subroutine Metadata__get_string

subroutine Metadata__set_array_int32(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: value(:)
  call atlas__Metadata__set_array_int(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_int32

subroutine Metadata__set_array_int64(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: value(:)
  call atlas__Metadata__set_array_long(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_int64

subroutine Metadata__set_array_real32(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value(:)
  call atlas__Metadata__set_array_float(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_real32

subroutine Metadata__set_array_real64(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value(:)
  call atlas__Metadata__set_array_double(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_real64

subroutine Metadata__get_array_int32(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_int), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_int(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Metadata__get_array_int32

subroutine Metadata__get_array_int64(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_long), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_long(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Metadata__get_array_int64

subroutine Metadata__get_array_real32(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_float), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_float(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Metadata__get_array_real32

subroutine Metadata__get_array_real64(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_double), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_double(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Metadata__get_array_real64


subroutine Metadata__get_mesh(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Mesh), intent(out) :: value
  value%cpp_object_ptr = atlas__Metadata__get_mesh(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_mesh

subroutine Metadata__get_grid(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_ReducedGrid), intent(out) :: value
  value%cpp_object_ptr = atlas__Metadata__get_grid(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_grid

subroutine Metadata__get_function_space(this, name, value)
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_FunctionSpace), intent(out) :: value
  value%cpp_object_ptr = atlas__Metadata__get_function_space(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_function_space

subroutine MetaData__print(this,channel)
  class(atlas_Metadata), intent(in) :: this
  class(atlas_LogChannel), intent(inout) :: channel
  call atlas__Metadata__print(this%cpp_object_ptr,channel%cpp_object_ptr)
end subroutine Metadata__print

function Metadata__json(this) result(json)
  character(len=:), allocatable :: json
  class(atlas_Metadata), intent(in) :: this
  type(c_ptr) :: json_cptr
  integer(c_int) :: json_size
  integer(c_int) :: json_allocated
  call atlas__Metadata__json(this%cpp_object_ptr,json_cptr,json_size,json_allocated)
  allocate(character(len=json_size) :: json )
  json = c_to_f_string_cptr(json_cptr)
  if( json_allocated == 1 ) call atlas_free(json_cptr)
end function Metadata__json

