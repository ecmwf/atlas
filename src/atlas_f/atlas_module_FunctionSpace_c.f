! (C) Copyright 2013-2015 ECMWF.

! -----------------------------------------------------------------------------
! FunctionSpace routines

function FunctionSpace__metadata(this) result(metadata)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  call metadata%reset_c_ptr( atlas__FunctionSpace__metadata(this%c_ptr()) )
end function FunctionSpace__metadata

subroutine FunctionSpace__create_field(this,name,nvars,kind)
  class(atlas_FunctionSpace), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: nvars
  integer, intent(in), optional :: kind

  ! ----------------------------------------------------
  ! ENUM DataType
  integer, parameter :: ATLAS_KIND_INT32  = -4
  integer, parameter :: ATLAS_KIND_INT64  = -8
  integer, parameter :: ATLAS_KIND_REAL32 =  4
  integer, parameter :: ATLAS_KIND_REAL64 =  8
  ! ----------------------------------------------------


  if (present(kind)) then
    if (kind == ATLAS_KIND_REAL64) then
      call atlas__FunctionSpace__create_field_double(this%c_ptr(),c_str(name),nvars)
    else if (kind == ATLAS_KIND_REAL32) then
      call atlas__FunctionSpace__create_field_float(this%c_ptr(),c_str(name),nvars)
    else if (kind == ATLAS_KIND_INT32) then
      call atlas__FunctionSpace__create_field_int(this%c_ptr(),c_str(name),nvars)
    else if (kind == ATLAS_KIND_INT64) then
      call atlas__FunctionSpace__create_field_long(this%c_ptr(),c_str(name),nvars)
    else
      write(0,*) "Unsupported kind"
      write(0,*) 'call abort()'
    endif
  else if (wp == c_double) then
    call atlas__FunctionSpace__create_field_double(this%c_ptr(),c_str(name),nvars)
  else if (wp == c_float) then
    call atlas__FunctionSpace__create_field_float(this%c_ptr(),c_str(name),nvars)
  end if
end subroutine FunctionSpace__create_field


subroutine FunctionSpace__remove_field(this,name)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__FunctionSpace__remove_field(this%c_ptr(),c_str(name))
end subroutine FunctionSpace__remove_field

function FunctionSpace__name(this) result(name)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__FunctionSpace__name(this%c_ptr())
  name = c_to_f_string_cptr(name_c_str)
end function FunctionSpace__name

function FunctionSpace__dof(this) result(dof)
  class(atlas_FunctionSpace), intent(in) :: this
  integer :: dof
  dof = atlas__FunctionSpace__dof(this%c_ptr())
end function FunctionSpace__dof

function FunctionSpace__glb_dof(this) result(glb_dof)
  class(atlas_FunctionSpace), intent(in) :: this
  integer :: glb_dof
  glb_dof = atlas__FunctionSpace__glb_dof(this%c_ptr())
end function FunctionSpace__glb_dof

function FunctionSpace__shape_arr(this) result(shape)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, pointer :: shape(:)
  type(c_ptr) :: shape_c_ptr
  integer(c_int) :: field_rank
  call atlas__FunctionSpace__shapef(this%c_ptr(), shape_c_ptr, field_rank)
  call C_F_POINTER ( shape_c_ptr , shape , (/field_rank/) )
end function FunctionSpace__shape_arr

function FunctionSpace__shape_idx(this,idx) result(shape_idx)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(in) :: idx
  integer :: shape_idx
  integer, pointer :: shape(:)
  shape => this%shape_arr()
  shape_idx = shape(idx)
end function FunctionSpace__shape_idx

function FunctionSpace__field(this,name) result(field)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FunctionSpace__field(this%c_ptr(), c_str(name) ) )
  if( field%is_null() ) write(0,*) 'call abort()'
  call field%return()
end function FunctionSpace__field

function FunctionSpace__has_field(this,name) result(flag)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FunctionSpace__has_field(this%c_ptr(), c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FunctionSpace__has_field

subroutine FunctionSpace__parallelise(this)
  class(atlas_FunctionSpace), intent(in) :: this
  call atlas__FunctionSpace__parallelise(this%c_ptr())
end subroutine FunctionSpace__parallelise

function FunctionSpace__get_halo_exchange(this) result(halo_exchange)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_HaloExchange) :: halo_exchange
  call halo_exchange%reset_c_ptr( atlas__FunctionSpace__halo_exchange( this%c_ptr() ) )
end function FunctionSpace__get_halo_exchange

subroutine FunctionSpace__halo_exchange_int32_r1(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_int( this%c_ptr(), field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r1
subroutine FunctionSpace__halo_exchange_int32_r2(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r2
subroutine FunctionSpace__halo_exchange_int32_r3(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r3

subroutine FunctionSpace__halo_exchange_real32_r1(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_float( this%c_ptr(), field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r1
subroutine FunctionSpace__halo_exchange_real32_r2(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r2
subroutine FunctionSpace__halo_exchange_real32_r3(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r3

subroutine FunctionSpace__halo_exchange_real64_r1(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_double( this%c_ptr(), field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r1
subroutine FunctionSpace__halo_exchange_real64_r2(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r2
subroutine FunctionSpace__halo_exchange_real64_r3(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r3
subroutine FunctionSpace__halo_exchange_real64_r4(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r4


function FunctionSpace__get_gather(this) result(gather)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_GatherScatter) :: gather
  call gather%reset_c_ptr( atlas__FunctionSpace__gather( this%c_ptr() ) )
end function FunctionSpace__get_gather

subroutine FunctionSpace__gather_real32_r1(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(in) :: field_data(:)
  real(c_float), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_float( this%c_ptr(), field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r1


subroutine FunctionSpace__gather_real32_r2(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_float( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r2


subroutine FunctionSpace__gather_real32_r3(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_float( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r3

subroutine FunctionSpace__gather_real64_r1(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(in) :: field_data(:)
  real(c_double), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_double( this%c_ptr(), field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r1


subroutine FunctionSpace__gather_real64_r2(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_double( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r2


subroutine FunctionSpace__gather_real64_r3(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_double( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r3

subroutine FunctionSpace__gather_int32_r1(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(in) :: field_data(:)
  integer, intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_int( this%c_ptr(), field_data, size(field_data), &
                                       & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r1

subroutine FunctionSpace__gather_int32_r2(this, field_data, glbfield_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(in) :: field_data(:,:)
  integer, intent(inout) :: glbfield_data(:,:)
  integer, pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_int( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r2

function FunctionSpace__get_checksum(this) result(checksum)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_Checksum) :: checksum
  call checksum%reset_c_ptr( atlas__FunctionSpace__checksum( this%c_ptr() ) )
end function FunctionSpace__get_checksum

subroutine atlas_FunctionSpace__delete(this)
  class(atlas_FunctionSpace), intent(inout) :: this
end subroutine


subroutine atlas_FunctionSpace__copy(this,obj_in)
  class(atlas_FunctionSpace), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine
