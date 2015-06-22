! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! Parametrisation routines

function atlas_Parametrisation__ctor() result(Parametrisation)
  use atlas_parametrisation_c_binding
  type(atlas_Parametrisation) :: Parametrisation
  Parametrisation%cpp_object_ptr = atlas__Parametrisation__new()
end function atlas_Parametrisation__ctor

subroutine atlas_Parametrisation__delete(this)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Parametrisation__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = C_NULL_ptr
end subroutine atlas_Parametrisation__delete

function Parametrisation__has(this, name) result(value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical :: value
  integer :: value_int
  value_int =  atlas__Parametrisation__has(this%cpp_object_ptr, c_str(name) )
  if( value_int == 1 ) then
    value = .True.
  else
    value = .False.
  end if
end function Parametrisation__has

subroutine Parametrisation__set_logical(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical, intent(in) :: value
  integer :: value_int
  if( value ) then
    value_int = 1
  else
    value_int = 0
  end if
  call atlas__Parametrisation__set_int(this%cpp_object_ptr, c_str(name), value_int )
end subroutine Parametrisation__set_logical

subroutine Parametrisation__set_int32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call atlas__Parametrisation__set_int(this%cpp_object_ptr, c_str(name), value)
end subroutine Parametrisation__set_int32

subroutine Parametrisation__set_real32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Parametrisation__set_float(this%cpp_object_ptr, c_str(name) ,value)
end subroutine Parametrisation__set_real32

subroutine Parametrisation__set_real64(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Parametrisation__set_double(this%cpp_object_ptr, c_str(name) ,value)
end subroutine Parametrisation__set_real64

subroutine Parametrisation__set_string(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Parametrisation__set_string(this%cpp_object_ptr, c_str(name) , c_str(value) )
end subroutine Parametrisation__set_string

subroutine Parametrisation__get_logical(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(out) :: value
  integer :: value_int
  value_int = atlas__Parametrisation__get_int(this%cpp_object_ptr,c_str(name) )
  if (value_int > 0) then
    value = .True.
  else
    value = .False.
  end if
end subroutine Parametrisation__get_logical

subroutine Parametrisation__get_int32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(out) :: value
  value = atlas__Parametrisation__get_int(this%cpp_object_ptr, c_str(name) )
end subroutine Parametrisation__get_int32

subroutine Parametrisation__get_real32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(out) :: value
  value = atlas__Parametrisation__get_float(this%cpp_object_ptr, c_str(name) )
end subroutine Parametrisation__get_real32

subroutine Parametrisation__get_real64(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(out) :: value
  value = atlas__Parametrisation__get_double(this%cpp_object_ptr, c_str(name) )
end subroutine Parametrisation__get_real64

subroutine Parametrisation__get_string(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(out) :: value
  character(len=MAX_STR_LEN) :: value_cstr
  call atlas__Parametrisation__get_string(this%cpp_object_ptr, c_str(name), value_cstr, MAX_STR_LEN )
  value = c_to_f_string_str(value_cstr)
end subroutine Parametrisation__get_string

subroutine Parametrisation__set_array_int32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: value(:)
  call atlas__Parametrisation__set_array_int(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Parametrisation__set_array_int32

subroutine Parametrisation__set_array_int64(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: value(:)
  call atlas__Parametrisation__set_array_long(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Parametrisation__set_array_int64

subroutine Parametrisation__set_array_real32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value(:)
  call atlas__Parametrisation__set_array_float(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Parametrisation__set_array_real32

subroutine Parametrisation__set_array_real64(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value(:)
  call atlas__Parametrisation__set_array_double(this%cpp_object_ptr, c_str(name), &
    & value, size(value) )
end subroutine Parametrisation__set_array_real64

subroutine Parametrisation__get_array_int32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_int), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Parametrisation__get_array_int(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Parametrisation__get_array_int32

subroutine Parametrisation__get_array_int64(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_long), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Parametrisation__get_array_long(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Parametrisation__get_array_int64

subroutine Parametrisation__get_array_real32(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_float), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Parametrisation__get_array_float(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Parametrisation__get_array_real32

subroutine Parametrisation__get_array_real64(this, name, value)
  use atlas_parametrisation_c_binding
  class(atlas_Parametrisation), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_double), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Parametrisation__get_array_double(this%cpp_object_ptr, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call atlas_free(value_cptr)
end subroutine Parametrisation__get_array_real64

