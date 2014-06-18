! (C) Copyright 2013-2014 ECMWF.


function new_Mesh() result(mesh)
  type(Mesh_type) :: mesh
  mesh%private%object = atlas__Mesh__new()
end function new_Mesh

subroutine Mesh__add_function_space(this, function_space)
  class(Mesh_type), intent(inout) :: this
  type(FunctionSpace_type), intent(in) :: function_space
  call atlas__Mesh__add_function_space(this%private%object,function_space%private%object)
end subroutine Mesh__add_function_space

function Mesh__function_space(this,name) result(function_space)
  class(Mesh_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(FunctionSpace_type) :: function_space
  function_space%private%object = atlas__Mesh__function_space(this%private%object, c_str(name) )
  if( .not. C_associated(function_space%private%object) ) write(0,*) 'call abort()'
end function Mesh__function_space

subroutine Mesh__delete(this)
  type(Mesh_type), intent(inout) :: this
  if ( c_associated(this%private%object) ) then
    call atlas__Mesh__delete(this%private%object)
  end if
  this%private%object = c_null_ptr
end subroutine Mesh__delete
