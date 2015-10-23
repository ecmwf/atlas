! (C) Copyright 2013-2015 ECMWF.

function atlas_Nabla__cptr(cptr) result(nabla)
  type(atlas_Nabla) :: nabla
  type(c_ptr), intent(in) :: cptr
  call nabla%reset_c_ptr( cptr )
end function

function atlas_Nabla__functionspace_config(functionspace,config) result(nabla)
  use atlas_nabla_c_binding
  type(atlas_Nabla) :: nabla
  class(atlas_NextFunctionSpace), intent(in) :: functionspace
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    nabla = atlas_Nabla__cptr(atlas__Nabla__create(functionspace%c_ptr(),config%c_ptr()))
  else
    opt_config = atlas_Config()
    nabla = atlas_Nabla__cptr(atlas__Nabla__create(functionspace%c_ptr(),opt_config%c_ptr()))
    call opt_config%final()
  endif
  call nabla%return()
end function

subroutine atlas_Nabla__delete(this)
  use atlas_nabla_c_binding
  class(atlas_Nabla), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Nabla__delete(this%c_ptr())
  endif
  call this%reset_c_ptr()
end subroutine atlas_Nabla__delete

subroutine atlas_Nabla__gradient(this,scalar,grad)
  use atlas_nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: scalar
  class(atlas_Field), intent(inout) :: grad
  call atlas__Nabla__gradient(this%c_ptr(),scalar%c_ptr(),grad%c_ptr())
end subroutine

subroutine atlas_Nabla__divergence(this,vector,div)
  use atlas_nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: vector
  class(atlas_Field), intent(inout) :: div
  call atlas__Nabla__divergence(this%c_ptr(),vector%c_ptr(),div%c_ptr())
end subroutine

subroutine atlas_Nabla__curl(this,vector,curl)
  use atlas_nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: vector
  class(atlas_Field), intent(inout) :: curl
  call atlas__Nabla__curl(this%c_ptr(),vector%c_ptr(),curl%c_ptr())
end subroutine

subroutine atlas_Nabla__laplacian(this,scalar,lapl)
  use atlas_nabla_c_binding
  class(atlas_Nabla), intent(in) :: this
  class(atlas_Field), intent(in) :: scalar
  class(atlas_Field), intent(inout) :: lapl
  call atlas__Nabla__laplacian(this%c_ptr(),scalar%c_ptr(),lapl%c_ptr())
end subroutine

! -----------------------------------------------------------------------------

