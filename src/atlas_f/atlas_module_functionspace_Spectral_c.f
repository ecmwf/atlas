! (C) Copyright 2013-2015 ECMWF.

! -----------------------------------------------------------------------------------------------

function atlas_functionspace_Spectral__cptr(cptr) result(functionspace)
  type(atlas_functionspace_Spectral) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function


#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_functionspace_Spectral__final(this)
  type(atlas_functionspace_Spectral), intent(inout) :: this
  call this%final()
end subroutine
#endif

function atlas_functionspace_Spectral__truncation(truncation) result(functionspace)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: functionspace
  integer(c_int), intent(in) :: truncation
  functionspace = atlas_functionspace_Spectral__cptr( &
    & atlas__SpectralFunctionSpace__new__truncation(truncation) )
  call functionspace%return()
end function

function atlas_functionspace_Spectral__trans(trans) result(functionspace)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: functionspace
  type(atlas_Trans), intent(in) :: trans
  functionspace = atlas_functionspace_Spectral__cptr( &
    & atlas__SpectralFunctionSpace__new__trans(trans%c_ptr()) )
  call functionspace%return()
end function

function SpectralFunctionSpace__create_field_name(this,name) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__SpectralFunctionSpace__create_field(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function SpectralFunctionSpace__create_field_name_lev(this,name,levels) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__SpectralFunctionSpace__create_field_lev(this%c_ptr(),c_str(name),levels) )
  call field%return()
end function

function SpectralFunctionSpace__create_glb_field_name(this,name) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__SpectralFunctionSpace__create_global_field(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function SpectralFunctionSpace__create_glb_field_name_lev(this,name,levels) result(field)
  use atlas_functionspace_spectral_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Spectral), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__SpectralFunctionSpace__create_global_field_lev(this%c_ptr(),c_str(name),levels) )
  call field%return()
end function

subroutine SpectralFunctionspace__gather(this,local,global)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__SpectralFunctionSpace__gather(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine SpectralFunctionspace__scatter(this,global,local)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__SpectralFunctionSpace__scatter(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

