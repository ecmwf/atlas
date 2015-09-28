! (C) Copyright 2013-2015 ECMWF.


function atlas_Mesh__ctor() result(mesh)
  type(atlas_Mesh) :: mesh
  mesh%cpp_object_ptr = atlas__Mesh__new()
end function atlas_Mesh__ctor

subroutine Mesh__create_function_space_nodes(this,name,shape_func,nb_nodes)
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_nodes
  integer :: shape(2)
  integer , parameter :: fortran_ordering = 1
  shape = (/FIELD_NB_VARS,nb_nodes/)
  call atlas__Mesh__create_function_space(this%cpp_object_ptr,c_str(name),c_str(shape_func), &
  & shape, size(shape), fortran_ordering )
end subroutine Mesh__create_function_space_nodes

subroutine Mesh__create_function_space_shape(this,name,shape_func,shape)
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: shape(:)
  integer , parameter :: fortran_ordering = 1
  call atlas__Mesh__create_function_space(this%cpp_object_ptr,c_str(name),c_str(shape_func), &
  & shape, size(shape), fortran_ordering )
end subroutine Mesh__create_function_space_shape

function Mesh__function_space(this,name) result(function_space)
  class(atlas_Mesh), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_FunctionSpace) :: function_space
  function_space%cpp_object_ptr = atlas__Mesh__function_space(this%cpp_object_ptr, c_str(name) )
  if( .not. C_associated(function_space%cpp_object_ptr) ) write(0,*) 'call abort()'
end function Mesh__function_space

function Mesh__create_nodes(this,nb_nodes) result(nodes)
  type(atlas_Nodes) :: nodes
  class(atlas_Mesh), intent(in) :: this
  integer, intent(in) :: nb_nodes
  nodes%cpp_object_ptr = atlas__Mesh__create_nodes(this%cpp_object_ptr,nb_nodes)
  if( .not. C_associated(nodes%cpp_object_ptr) ) write(0,*) 'call abort()'
end function

function Mesh__nodes(this) result(nodes)
  class(atlas_Mesh), intent(in) :: this
  type(atlas_Nodes) :: nodes
  nodes%cpp_object_ptr = atlas__Mesh__nodes(this%cpp_object_ptr)
  if( .not. C_associated(nodes%cpp_object_ptr) ) write(0,*) 'call abort()'
end function

subroutine atlas_Mesh__delete(this)
  type(atlas_Mesh), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Mesh__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine atlas_Mesh__delete

subroutine atlas_Mesh_finalize( this )
  class(atlas_Mesh), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__Mesh__delete(this%cpp_object_ptr);
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine

