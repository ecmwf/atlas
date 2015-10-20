! (C) Copyright 2013-2015 ECMWF.

function atlas_functionspace_EdgeBasedFiniteVolume__cptr(cptr) result(functionspace)
  type(atlas_functionspace_EdgeBasedFiniteVolume) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function

function atlas_functionspace_EdgeBasedFiniteVolume__mesh_halo(mesh,halo) result(function_space)
  use atlas_functionspace_EdgeBasedFiniteVolume_c_binding
  type(atlas_functionspace_EdgeBasedFiniteVolume) :: function_space
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  integer :: opt_halo
  opt_halo = 1
  if( present(halo) ) opt_halo = halo
  function_space = atlas_functionspace_EdgeBasedFiniteVolume__cptr( &
    & atlas__functionspace__EdgeBasedFiniteVolume__new(mesh%c_ptr(),opt_halo) )
  call function_space%return()
end function

! -----------------------------------------------------------------------------------------------

