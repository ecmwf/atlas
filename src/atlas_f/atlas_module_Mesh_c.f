! (C) Copyright 2013-2015 ECMWF.

function atlas_Mesh__cptr(cptr) result(mesh)
  type(atlas_Mesh) :: mesh
  type(c_ptr), intent(in) :: cptr
  call mesh%reset_c_ptr( cptr )
end function atlas_Mesh__cptr

function atlas_Mesh__ctor() result(mesh)
  type(atlas_Mesh) :: mesh
  call mesh%reset_c_ptr( atlas__Mesh__new() )
end function atlas_Mesh__ctor

subroutine Mesh__create_function_space_nodes(this,name,shape_func,nb_nodes)
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_nodes
  integer :: shape(2)
  integer , parameter :: fortran_ordering = 1
  shape = (/FIELD_NB_VARS,nb_nodes/)
  call atlas__Mesh__create_function_space(this%c_ptr(),c_str(name),c_str(shape_func), &
  & shape, size(shape), fortran_ordering )
end subroutine Mesh__create_function_space_nodes

subroutine Mesh__create_function_space_shape(this,name,shape_func,shape)
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: shape(:)
  integer , parameter :: fortran_ordering = 1
  call atlas__Mesh__create_function_space(this%c_ptr(),c_str(name),c_str(shape_func), &
  & shape, size(shape), fortran_ordering )
end subroutine Mesh__create_function_space_shape

function Mesh__function_space(this,name) result(function_space)
  class(atlas_Mesh), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_deprecated_FunctionSpace) :: function_space
  call function_space%reset_c_ptr( atlas__Mesh__function_space(this%c_ptr(), c_str(name) ) )
  if( function_space%is_null() ) write(0,*) 'call abort()'
end function Mesh__function_space

function Mesh__create_nodes(this,nb_nodes) result(nodes)
  type(atlas_mesh_Nodes) :: nodes
  class(atlas_Mesh), intent(in) :: this
  integer, intent(in) :: nb_nodes
  call nodes%reset_c_ptr( atlas__Mesh__create_nodes(this%c_ptr(),nb_nodes) )
  if( nodes%is_null() ) write(0,*) 'call abort()'
end function

function Mesh__nodes(this) result(nodes)
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Nodes) :: nodes
  call nodes%reset_c_ptr( atlas__Mesh__nodes(this%c_ptr()) )
  if( nodes%is_null() ) write(0,*) 'call abort()'
end function

function Mesh__cells(this) result(cells)
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Cells) :: cells
  call cells%reset_c_ptr( atlas__Mesh__cells(this%c_ptr()) )
  if( cells%is_null() ) write(0,*) 'call abort()'
end function

function Mesh__edges(this) result(cells)
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Edges) :: cells
  call cells%reset_c_ptr( atlas__Mesh__edges(this%c_ptr()) )
  if( cells%is_null() ) write(0,*) 'call abort()'
end function

subroutine atlas_Mesh__delete(this)
  class(atlas_Mesh), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Mesh__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_Mesh__delete

subroutine atlas_Mesh__copy(this,obj_in)
  class(atlas_Mesh), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine
