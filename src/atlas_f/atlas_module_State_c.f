! (C) Copyright 2013-201 ECMWF.

function atlas_State__new() result(State)
  use atlas_state_c_binding
  type(atlas_State) :: State
  State%cpp_object_ptr = atlas__State__new()
end function

function atlas_State__create(state_type, params) result(State)
  use atlas_state_c_binding
  type(atlas_State) :: State
  character(len=*), intent(in) :: state_type
  class(atlas_Parametrisation), intent(in) :: params
  State%cpp_object_ptr = atlas__State__create(c_str(state_type),params%cpp_object_ptr)
end function

subroutine atlas_State__delete(this)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__State__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine

subroutine atlas_State__add_field(this,field)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  class(atlas_Field), intent(in) :: field
  call atlas__State__add_field(this%cpp_object_ptr,field%cpp_object_ptr)
end subroutine

subroutine atlas_State__remove_field(this,name)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  call atlas__State__remove_field(this%cpp_object_ptr,c_str(name))
end subroutine

function atlas_State__has_field(this,name) result(has_field)
  use atlas_state_c_binding
  logical :: has_field
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer :: has_field_int
  has_field_int = atlas__State__has_field(this%cpp_object_ptr,c_str(name))
  has_field = .False.
  if( has_field_int == 1 ) has_field = .True.
end function

function atlas_State__nb_fields(this) result(nb_fields)
  use atlas_state_c_binding
  integer :: nb_fields
  class(atlas_State), intent(inout) :: this
  nb_fields = atlas__State__nb_fields(this%cpp_object_ptr)
end function

function atlas_State__field_by_name(this,name) result(field)
  use atlas_state_c_binding
  type(atlas_Field) :: field
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  field%cpp_object_ptr = atlas__State__field_by_name(this%cpp_object_ptr,c_str(name))
end function

function atlas_State__field_by_index(this,index) result(field)
  use atlas_state_c_binding
  type(atlas_Field) :: field
  class(atlas_State), intent(inout) :: this
  integer, intent(in) :: index
  field%cpp_object_ptr = atlas__State__field_by_index(this%cpp_object_ptr,index-1)
end function
