! (C) Copyright 2013-2015 ECMWF.

! -----------------------------------------------------------------------------------------------

function atlas_SpectralFunctionSpace__cptr(cptr) result(functionspace)
  type(atlas_SpectralFunctionSpace) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
  call functionspace%attach()
  call atlas_return(functionspace)
end function

subroutine atlas_SpectralFunctionSpace__finalize(this)
  use atlas_Spectralfunctionspace_c_binding
  class(atlas_SpectralFunctionSpace), intent(inout) :: this
  if( .not. this%is_null() ) then
    if( this%owners() <= 0 ) then
      call atlas_abort("Cannot finalize functionspace that has no owners")
    endif
    call this%detach()
    if( this%owners() == 0 ) then
      call atlas__SpectralFunctionSpace__delete(this%c_ptr())
    endif
    call this%reset_c_ptr()
  endif
end subroutine

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_SpectralFunctionSpace__final(this)
  type(atlas_SpectralFunctionSpace), intent(inout) :: this
  call this%finalize()
end subroutine
#endif

subroutine atlas_SpectralFunctionSpace__delete(this)
  use atlas_Spectralfunctionspace_c_binding
  type(atlas_SpectralFunctionSpace), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__SpectralFunctionSpace__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_SpectralFunctionSpace__delete

subroutine atlas_SpectralFunctionSpace__reset(functionspace_out,functionspace_in)
  type(atlas_SpectralFunctionSpace), intent(inout) :: functionspace_out
  class(atlas_SpectralFunctionSpace), intent(in) :: functionspace_in
  if( .not. atlas_compare_equal(functionspace_out%c_ptr(),functionspace_in%c_ptr()) ) then
#ifndef FORTRAN_SUPPORTS_FINAL
    call atlas_SpectralFunctionSpace__finalize(functionspace_out)
#endif
    call functionspace_out%reset_c_ptr( functionspace_in%c_ptr() )
    if( .not. functionspace_out%is_null() ) call functionspace_out%attach()
  endif
end subroutine


function atlas_SpectralFunctionSpace__name_truncation(name,truncation) result(functionspace)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_SpectralFunctionSpace) :: functionspace
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: truncation
  functionspace = atlas_SpectralFunctionSpace__cptr( &
    & atlas__SpectralFunctionSpace__new__name_truncation(c_str(name),truncation) )
  call atlas_return(functionspace)
end function

function atlas_SpectralFunctionSpace__name_trans(name,trans) result(functionspace)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_SpectralFunctionSpace) :: functionspace
  character(len=*), intent(in) :: name
  type(atlas_Trans), intent(in) :: trans
  functionspace = atlas_SpectralFunctionSpace__cptr( &
    & atlas__SpectralFunctionSpace__new__name_trans(c_str(name),trans%c_ptr()) )
  call atlas_return(functionspace)
end function

function SpectralFunctionSpace__create_field_name(this,name) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__SpectralFunctionSpace__create_field(this%c_ptr(),c_str(name)) )
  call atlas_return(field)
end function

function SpectralFunctionSpace__create_field_name_lev(this,name,levels) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__SpectralFunctionSpace__create_field_lev(this%c_ptr(),c_str(name),levels) )
  call atlas_return(field)
end function

function SpectralFunctionSpace__create_glb_field_name(this,name) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__SpectralFunctionSpace__create_global_field(this%c_ptr(),c_str(name)) )
  call atlas_return(field)
end function

function SpectralFunctionSpace__create_glb_field_name_lev(this,name,levels) result(field)
  use atlas_spectralfunctionspace_c_binding
  type(atlas_Field) :: field
  class(atlas_SpectralFunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__SpectralFunctionSpace__create_global_field_lev(this%c_ptr(),c_str(name),levels) )
  call atlas_return(field)
end function

