! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! Field routines

function Field__name(this) result(field_name)
  class(Field_type), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = atlas__Field__name(this%private%object)
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

function Field__data_type(this) result(field_data_type)
  class(Field_type), intent(in) :: this
  character(len=:), allocatable :: field_data_type
  type(c_ptr) :: field_data_type_c_str
  field_data_type_c_str = atlas__Field__data_type(this%private%object)
  field_data_type = c_to_f_string_cptr(field_data_type_c_str)
end function Field__data_type

function Field__nb_vars(this) result(nb_vars)
  class(Field_type), intent(in) :: this
  integer :: nb_vars
  nb_vars = atlas__Field__nb_vars(this%private%object)
end function Field__nb_vars

function Field__metadata(this) result(metadata)
  class(Field_type), intent(in) :: this
  type(metadata_type) :: Metadata
  metadata%private%object = atlas__Field__metadata(this%private%object)
end function Field__metadata

function Field__function_space(this) result(function_space)
  class(Field_type), intent(in) :: this
  type(FunctionSpace_type) :: function_space
  function_space%private%object = atlas__Field__function_space(this%private%object)
end function Field__function_space

subroutine Field__access_data1_integer(this, field) 
  class(Field_type), intent(in) :: this
  integer, pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_boundsf_int(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_integer

subroutine Field__access_data2_integer(this, field) 
  class(Field_type), intent(in) :: this
  integer, pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_int(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_integer

subroutine Field__access_data3_integer(this, field) 
  class(Field_type), intent(in) :: this
  integer, pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_int(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_integer

subroutine Field__access_data1_real32(this, field) 
  class(Field_type), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_boundsf_float(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_real32

subroutine Field__access_data2_real32(this, field) 
  class(Field_type), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_float(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_real32

subroutine Field__access_data3_real32(this, field) 
  class(Field_type), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_float(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_real32

subroutine Field__access_data1_real64(this, field) 
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_real64

subroutine Field__access_data2_real64(this, field) 
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_real64

subroutine Field__access_data3_real64(this, field) 
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_real64

subroutine Field__access_data3_real64_bounds(this, field, field_bounds)
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:)
  integer, intent(in) :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end subroutine Field__access_data3_real64_bounds

function Field__data1_wp(this) result(field)
  class(Field_type), intent(in) :: this
  real(wp), pointer :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  if( wp == c_double ) then
    call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  else if (wp == c_float ) then
    call atlas__Field__data_boundsf_float(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  end if
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end function Field__data1_wp

function Field__data2_wp(this) result(field)
  class(Field_type), intent(in) :: this
  real(wp), pointer :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  if( wp == c_double ) then
    call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  else if (wp == c_float ) then
    call atlas__Field__data_boundsf_float(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  end if
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  if (size(field_bounds) < 2) then
    write(0,*) "Cannot access field """,this%name(),""" with rank",field_rank," as rank 2"
    write(0,*) 'call abort()'
  end if
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
  if( size(field) /= field_size ) then
    write(0,*) "Requested bounds of field ", this%name(), "[", field_bounds(1:2), &
     & "] do not cover the entire field of size ", field_size
    write(0,*) 'call abort()'
  end if
end function Field__data2_wp

function Field__data3_wp(this) result(field)
  class(Field_type), intent(in) :: this
  real(wp), pointer :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  if( wp == c_double ) then
    call atlas__Field__data_boundsf_double(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  else if (wp == c_float ) then
    call atlas__Field__data_boundsf_float(this%private%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  end if
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  if (size(field_bounds) < 3) then
    write(0,*) "Cannot access field """,this%name(),""" with rank",field_rank," as rank 3"
    write(0,*) 'call abort()'
  end if
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , field_bounds(1:3) )
  if( size(field) /= field_size ) then
    write(0,*) "Requested bounds of field ", field_bounds(1:3), " do not cover the entire field of size ", field_size
    write(0,*) 'call abort()'
  end if
end function Field__data3_wp
