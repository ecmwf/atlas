! (C) Copyright 2013-201 ECMWF.

function atlas_State__new() result(State)
  use atlas_state_c_binding
  type(atlas_State) :: State
  State%cpp_object_ptr = atlas__State__new()
end function

function atlas_State__generate(generator, params) result(State)
  use atlas_state_c_binding
  type(atlas_State) :: State
  character(len=*), intent(in) :: generator
  class(atlas_Config), intent(in), optional :: params

  type(atlas_Config) :: p

  State%cpp_object_ptr = atlas__State__new()

  if( present(params) ) then
    call atlas__State__initialize(State%cpp_object_ptr,c_str(generator),params%cpp_object_ptr)
  else
    p = atlas_Config()
    call atlas__State__initialize(State%cpp_object_ptr,c_str(generator),p%cpp_object_ptr)
    call atlas_delete(p)
  endif
end function

subroutine atlas_State__delete(this)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__State__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine

subroutine atlas_State__add(this,field)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  class(atlas_Field), intent(in) :: field
  call atlas__State__add(this%cpp_object_ptr,field%cpp_object_ptr)
end subroutine

subroutine atlas_State__remove(this,name)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  call atlas__State__remove(this%cpp_object_ptr,c_str(name))
end subroutine

function atlas_State__has(this,name) result(has)
  use atlas_state_c_binding
  logical :: has
  class(atlas_State), intent(in) :: this
  character(len=*), intent(in) :: name
  integer :: has_int
  has_int = atlas__State__has(this%cpp_object_ptr,c_str(name))
  has = .False.
  if( has_int == 1 ) has = .True.
end function

function atlas_State__size(this) result(size)
  use atlas_state_c_binding
  integer :: size
  class(atlas_State), intent(in) :: this
  size = atlas__State__size(this%cpp_object_ptr)
end function

function atlas_State__field_by_name(this,name) result(field)
  use atlas_state_c_binding
  type(atlas_Field) :: field
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__State__field_by_name(this%cpp_object_ptr,c_str(name)) )
end function

function atlas_State__field_by_index(this,index) result(field)
  use atlas_state_c_binding
  type(atlas_Field) :: field
  class(atlas_State), intent(in) :: this
  integer, intent(in) :: index
  field = atlas_Field( atlas__State__field_by_index(this%cpp_object_ptr,index-1) )
end function

function atlas_State__metadata(this) result(metadata)
  use atlas_state_c_binding
  type(atlas_Metadata) :: metadata
  class(atlas_State), intent(in) :: this
  metadata%cpp_object_ptr = atlas__State__metadata(this%cpp_object_ptr)
end function


