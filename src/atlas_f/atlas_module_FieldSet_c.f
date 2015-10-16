! (C) Copyright 2013-2015 ECMWF.

! -----------------------------------------------------------------------------
! FieldSet routines

function atlas_FieldSet__ctor(name) result(fieldset)
  character(len=*), intent(in) :: name
  type(atlas_FieldSet) :: fieldset
  call fieldset%reset_c_ptr( atlas__FieldSet__new( c_str(name) ) )
end function atlas_FieldSet__ctor

subroutine atlas_FieldSet__delete(this)
  type(atlas_FieldSet), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__FieldSet__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_FieldSet__delete

subroutine FieldSet__add_field(this,field)
  class(atlas_FieldSet), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__FieldSet__add_field(this%c_ptr(), field%c_ptr())
end subroutine FieldSet__add_field

function FieldSet__has_field(this,name) result(flag)
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FieldSet__has_field(this%c_ptr(), c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FieldSet__has_field

function FieldSet__size(this) result(nb_fields)
  class(atlas_FieldSet), intent(in) :: this
  integer :: nb_fields
  nb_fields = atlas__FieldSet__size(this%c_ptr())
end function FieldSet__size

function FieldSet__field_by_name(this,name) result(field)
  class(atlas_FieldSet), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_name(this%c_ptr(), c_str(name) ) )
  call field%return()
end function FieldSet__field_by_name

function FieldSet__field_by_idx(this,idx) result(field)
  class(atlas_FieldSet), intent(in) :: this
  integer, intent(in) :: idx
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FieldSet__field_by_idx(this%c_ptr(), idx-1) ) ! C index
  call field%return()
end function FieldSet__field_by_idx
