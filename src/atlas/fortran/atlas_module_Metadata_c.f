! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! Metadata routines

function new_Metadata() result(metadata)
  type(Metadata_type) :: metadata
  metadata%cpp_object_ptr = atlas__Metadata__new()
end function new_Metadata

subroutine Metadata__delete(this)
  type(Metadata_type), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Metadata__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = C_NULL_ptr
end subroutine Metadata__delete

function Metadata__has(this, name) result(value)
  class(Metadata_type), intent(inout) :: this
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

subroutine Metadata__add_logical(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical, intent(in) :: value
  integer :: value_int
  if( value ) then
    value_int = 1
  else
    value_int = 0
  end if
  call atlas__Metadata__add_int(this%cpp_object_ptr, c_str(name), value_int )
end subroutine Metadata__add_logical

subroutine Metadata__add_integer(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call atlas__Metadata__add_int(this%cpp_object_ptr, c_str(name), value)
end subroutine Metadata__add_integer

subroutine Metadata__add_real32(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Metadata__add_float(this%cpp_object_ptr, c_str(name) ,value)
end subroutine Metadata__add_real32

subroutine Metadata__add_real64(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Metadata__add_double(this%cpp_object_ptr, c_str(name) ,value)
end subroutine Metadata__add_real64

subroutine Metadata__add_string(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Metadata__add_string(this%cpp_object_ptr, c_str(name) , c_str(value) )
end subroutine Metadata__add_string

subroutine Metadata__get_logical(this, name, value)
  class(Metadata_type), intent(in) :: this
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

subroutine Metadata__get_integer(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(out) :: value
  value = atlas__Metadata__get_int(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_integer

subroutine Metadata__get_real32(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(out) :: value
  value = atlas__Metadata__get_float(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_real32

subroutine Metadata__get_real64(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(out) :: value
  value = atlas__Metadata__get_double(this%cpp_object_ptr, c_str(name) )
end subroutine Metadata__get_real64

subroutine Metadata__get_string(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(out) :: value
  type(c_ptr) :: value_c_str
  value_c_str = atlas__Metadata__get_string(this%cpp_object_ptr, c_str(name) )
  value = c_to_f_string_cptr(value_c_str)
end subroutine Metadata__get_string
