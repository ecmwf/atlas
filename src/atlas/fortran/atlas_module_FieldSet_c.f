! (C) Copyright 2013-2014 ECMWF.

! -----------------------------------------------------------------------------
! FieldSet routines

function new_atlas_FieldSet(name) result(fieldset)
  character(len=*), intent(in) :: name
  type(atlas_FieldSet) :: fieldset
  fieldset%cpp_object_ptr = atlas__FieldSet__new( c_str(name) )
end function new_atlas_FieldSet

subroutine atlas_FieldSet__delete(this)
  type(atlas_FieldSet), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__FieldSet__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = C_NULL_ptr
end subroutine atlas_FieldSet__delete

subroutine FieldSet__add_field(this,field)
  class(atlas_FieldSet), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__FieldSet__add_field(this%cpp_object_ptr, field%cpp_object_ptr)
end subroutine FieldSet__add_field

function FieldSet__has_field(this,name) result(flag)
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FieldSet__has_field(this%cpp_object_ptr, c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FieldSet__has_field

function FieldSet__size(this) result(nb_fields)
  class(atlas_FieldSet), intent(in) :: this
  integer :: nb_fields
  nb_fields = atlas__FieldSet__size(this%cpp_object_ptr)
end function FieldSet__size

function FieldSet__field_by_name(this,name) result(field)
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field%cpp_object_ptr = atlas__FieldSet__field_by_name(this%cpp_object_ptr, c_str(name) )
end function FieldSet__field_by_name

function FieldSet__field_by_idx(this,idx) result(field)
  class(atlas_FieldSet), intent(in) :: this
  integer, intent(in) :: idx
  type(atlas_Field) :: field
  field%cpp_object_ptr = atlas__FieldSet__field_by_idx(this%cpp_object_ptr, idx-1) ! C index
end function FieldSet__field_by_idx

subroutine FieldSet__fields(this,fields)
  class(atlas_FieldSet), intent(in) :: this
  type(atlas_Field), allocatable, intent(out) :: fields(:)

  type(c_ptr), pointer :: fields_ptr(:)
  type(c_ptr) :: fields_cptr
  integer :: nb_fields, jfield
  fields_cptr = c_null_ptr
  call atlas__FieldSet__fields(this%cpp_object_ptr, fields_cptr, nb_fields)
  call c_f_pointer( fields_cptr, fields_ptr, (/nb_fields/) )
  allocate( fields(nb_fields) )
  do jfield=1,nb_fields
    fields(jfield)%cpp_object_ptr = fields_ptr(jfield)
  end do
end subroutine FieldSet__fields
