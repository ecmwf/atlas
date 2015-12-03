! (C) Copyright 2013-2015 ECMWF.

! -----------------------------------------------------------------------------------------------

function ReducedGridPoints__cptr(cptr) result(functionspace)
  type(atlas_functionspace_ReducedGridPoints) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function


#ifdef FORTRAN_SUPPORTS_FINAL
subroutine ReducedGridPoints__final(this)
  type(atlas_functionspace_ReducedGridPoints), intent(inout) :: this
  call this%final()
end subroutine
#endif

function ReducedGridPoints__grid(grid) result(functionspace)
  use atlas_functionspace_reducedgridpoints_c_binding
  type(atlas_functionspace_ReducedGridPoints) :: functionspace
  class(atlas_Grid), intent(in) :: grid
  functionspace = ReducedGridPoints__cptr( &
    & atlas__functionspace__ReducedGridPoints__new__grid( grid%c_ptr()) )
  call functionspace%return()
end function

function ReducedGridPoints__create_field_name(this,name) result(field)
  use atlas_functionspace_reducedgridpoints_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridPoints) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__functionspace__ReducedGridPoints__create_field(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function ReducedGridPoints__create_field_name_lev(this,name,levels) result(field)
  use atlas_functionspace_reducedgridpoints_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridPoints), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__functionspace__ReducedGridPoints__create_field_lev(this%c_ptr(),c_str(name),levels) )
  call field%return()
end function

function ReducedGridPoints__create_glb_field_name(this,name) result(field)
  use atlas_functionspace_reducedgridpoints_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridPoints) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__functionspace__ReducedGridPoints__create_gfield(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function ReducedGridPoints__create_glb_field_name_lev(this,name,levels) result(field)
  use atlas_functionspace_reducedgridpoints_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_ReducedGridPoints), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  field = atlas_Field( atlas__functionspace__ReducedGridPoints__create_gfield_lev(this%c_ptr(),c_str(name),levels) )
  call field%return()
end function

subroutine ReducedGridPoints__gather(this,local,global)
  use atlas_functionspace_reducedgridpoints_c_binding
  class(atlas_functionspace_ReducedGridPoints), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__functionspace__ReducedGridPoints__gather(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

subroutine ReducedGridPoints__scatter(this,global,local)
  use atlas_functionspace_reducedgridpoints_c_binding
  class(atlas_functionspace_ReducedGridPoints), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__functionspace__ReducedGridPoints__scatter(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

function ReducedGridPoints__checksum_fieldset(this,fieldset) result(checksum)
  use atlas_functionspace_reducedgridpoints_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_ReducedGridPoints), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__fs__ReducedGridPoints__checksum_fieldset(this%c_ptr(),fieldset%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


function ReducedGridPoints__checksum_field(this,field) result(checksum)
  use atlas_functionspace_reducedgridpoints_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_ReducedGridPoints), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__fs__ReducedGridPoints__checksum_field(this%c_ptr(),field%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function

