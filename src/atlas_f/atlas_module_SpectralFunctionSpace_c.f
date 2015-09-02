! (C) Copyright 2013-2015 ECMWF.

! -----------------------------------------------------------------------------------------------

function atlas_SpectralFunctionSpace__name_truncation(name,truncation) result(functionspace)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_SpectralFunctionSpace) :: functionspace
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: truncation
  functionspace%cpp_object_ptr = atlas__SpectralFunctionSpace__new__name_truncation(c_str(name),truncation)
end function

function atlas_SpectralFunctionSpace__name_trans(name,trans) result(functionspace)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_SpectralFunctionSpace) :: functionspace
  character(len=*), intent(in) :: name
  type(atlas_Trans), intent(in) :: trans
  functionspace%cpp_object_ptr = atlas__SpectralFunctionSpace__new__name_trans(c_str(name),trans%cpp_object_ptr)
end function

subroutine atlas_SpectralFunctionSpace__delete(this)
  use atlas_spectralfunctionspace_c_binding
  class(atlas_SpectralFunctionSpace) :: this
  call atlas__SpectralFunctionSpace__delete(This%cpp_object_ptr)
end subroutine

function SpectralFunctionSpace__create_field_name(this,name) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__SpectralFunctionSpace__create_field(this%cpp_object_ptr,c_str(name)) )
end function

function SpectralFunctionSpace__create_field_name_lev(this,name,levels) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__SpectralFunctionSpace__create_field_lev(this%cpp_object_ptr,c_str(name),levels) )
end function

function SpectralFunctionSpace__create_glb_field_name(this,name) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__SpectralFunctionSpace__create_global_field(this%cpp_object_ptr,c_str(name)) )
end function

function SpectralFunctionSpace__create_glb_field_name_lev(this,name,levels) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__SpectralFunctionSpace__create_global_field_lev(this%cpp_object_ptr,c_str(name),levels) )
end function

