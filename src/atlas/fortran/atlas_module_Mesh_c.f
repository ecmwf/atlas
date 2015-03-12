! (C) Copyright 2013-2014 ECMWF.


function new_atlas_Mesh() result(mesh)
  type(atlas_Mesh) :: mesh
  mesh%cpp_object_ptr = atlas__Mesh__new()
end function new_atlas_Mesh

subroutine Mesh__create_function_space(this,name,shape_func,nb_nodes)
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_nodes
  integer :: extents(2)
  extents = (/nb_nodes,FIELD_NB_VARS/)
  call atlas__Mesh__create_function_space(this%cpp_object_ptr,c_str(name),c_str(shape_func), &
  & extents, 2)
end subroutine Mesh__create_function_space

function Mesh__function_space(this,name) result(function_space)
  class(atlas_Mesh), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_FunctionSpace) :: function_space
  function_space%cpp_object_ptr = atlas__Mesh__function_space(this%cpp_object_ptr, c_str(name) )
  if( .not. C_associated(function_space%cpp_object_ptr) ) write(0,*) 'call abort()'
end function Mesh__function_space

subroutine atlas_Mesh__delete(this)
  type(atlas_Mesh), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Mesh__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine atlas_Mesh__delete
