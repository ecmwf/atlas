! (C) Copyright 2013-2014 ECMWF.


function atlas_NodesFunctionSpace__ctor(name,mesh,halo) result(function_space)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesFunctionSpace) :: function_space
  character(len=*), intent(in) :: name
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 0
  if( present(halo) ) opt_halo = halo
  function_space%cpp_object_ptr = atlas__NodesFunctionSpace__new(c_str(name),mesh%cpp_object_ptr,opt_halo)
end function atlas_NodesFunctionSpace__ctor

subroutine atlas_NodesFunctionSpace__delete(this)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesFunctionSpace), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__NodesFunctionSpace__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine atlas_NodesFunctionSpace__delete


function atlas_NodesColumnFunctionSpace__ctor(name,mesh,nb_levels,halo) result(function_space)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesColumnFunctionSpace) :: function_space
  character(len=*), intent(in) :: name
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in) :: nb_levels
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 0
  if( present(halo) ) opt_halo = halo
  function_space%cpp_object_ptr = atlas__NodesColumnFunctionSpace__new(c_str(name),mesh%cpp_object_ptr,nb_levels,opt_halo)
end function atlas_NodesColumnFunctionSpace__ctor

subroutine atlas_NodesColumnFunctionSpace__delete(this)
  use atlas_nodesfunctionspace_c_binding
  type(atlas_NodesColumnFunctionSpace), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__NodesColumnFunctionSpace__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine atlas_NodesColumnFunctionSpace__delete