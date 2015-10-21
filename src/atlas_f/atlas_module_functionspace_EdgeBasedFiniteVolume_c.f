! (C) Copyright 2013-2015 ECMWF.

function atlas_functionspace_EdgeBasedFiniteVolume__cptr(cptr) result(functionspace)
  type(atlas_functionspace_EdgeBasedFiniteVolume) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function

function atlas_functionspace_EdgeBasedFiniteVolume__mesh_config(mesh,config) result(function_space)
  use atlas_functionspace_EdgeBasedFiniteVolume_c_binding
  type(atlas_functionspace_EdgeBasedFiniteVolume) :: function_space
  type(atlas_Mesh), intent(inout) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    function_space = atlas_functionspace_EdgeBasedFiniteVolume__cptr( &
      & atlas__functionspace__EdgeBasedFiniteVolume__new(mesh%c_ptr(),config%c_ptr()) )
  else
    opt_config = atlas_Config()
    function_space = atlas_functionspace_EdgeBasedFiniteVolume__cptr( &
      & atlas__functionspace__EdgeBasedFiniteVolume__new(mesh%c_ptr(),opt_config%c_ptr()) )
    call opt_config%final()
  endif
  call function_space%return()
end function

! -----------------------------------------------------------------------------------------------

