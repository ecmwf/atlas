! (C) Copyright 2013-2015 ECMWF.

function atlas_MeshGenerator__cptr(cptr) result(MeshGenerator)
  type(atlas_MeshGenerator) :: MeshGenerator
  type(c_ptr), intent(in) :: cptr
  call MeshGenerator%reset_c_ptr( cptr )
end function

function atlas_MeshGenerator__name_config(name,config) result(MeshGenerator)
  use atlas_MeshGenerator_c_binding
  type(atlas_MeshGenerator) :: MeshGenerator
  character(len=*), intent(in) :: name
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    MeshGenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str(name),config%c_ptr()))
  else
    opt_config = atlas_Config()
    MeshGenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str(name),opt_config%c_ptr()))
    call opt_config%final()
  endif
  call MeshGenerator%return()
end function

subroutine atlas_MeshGenerator__delete(this)
  use atlas_MeshGenerator_c_binding
  class(atlas_MeshGenerator), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__MeshGenerator__delete(this%c_ptr())
  endif
  call this%reset_c_ptr()
end subroutine atlas_MeshGenerator__delete

function atlas_MeshGenerator__generate(this,grid,distribution) result(mesh)
   use atlas_MeshGenerator_c_binding
   type(atlas_Mesh) :: mesh
   class(atlas_MeshGenerator), intent(in) :: this
   class(atlas_ReducedGrid), intent(in) :: grid
   class(atlas_GridDistribution), intent(in), optional :: distribution

   if( present(distribution) ) then
     mesh = atlas_Mesh( atlas__MeshGenerator__generate__grid_griddist(this%c_ptr(),grid%c_ptr(),distribution%c_ptr()) )
   else
     mesh = atlas_Mesh( atlas__MeshGenerator__generate__grid(this%c_ptr(),grid%c_ptr()) )
   endif
   call mesh%return()
end function


! -----------------------------------------------------------------------------

