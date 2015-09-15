! (C) Copyright 2013-2015 ECMWF.

function atlas_NextFunctionSpace__cptr(cptr) result(functionspace)
  type(atlas_NextFunctionSpace) :: functionspace
  type(c_ptr), intent(in) :: cptr
  functionspace%cpp_object_ptr = cptr
  call functionspace%attach()
  call atlas_return(functionspace)
end function

subroutine atlas_NextFunctionSpace__finalize(this)
  use atlas_functionspace_c_binding
  class(atlas_NextFunctionSpace), intent(inout) :: this
  if( c_associated(this%cpp_object_ptr) ) then
    if( this%owners() <= 0 ) then
      call atlas_abort("Cannot finalize functionspace that has no owners")
    endif
    call this%detach()
    if( this%owners() == 0 ) then
      call atlas__NextFunctionSpace__delete(this%cpp_object_ptr)
    endif
    this%cpp_object_ptr = c_null_ptr
  endif
end subroutine

#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_NextFunctionSpace__final(this)
  type(atlas_NextFunctionSpace), intent(inout) :: this
  call this%finalize()
end subroutine
#endif

subroutine atlas_NextFunctionSpace__delete(this)
  use atlas_functionspace_c_binding
  type(atlas_NextFunctionSpace), intent(inout) :: this
  if ( c_associated(this%cpp_object_ptr) ) then
    call atlas__NextFunctionSpace__delete(this%cpp_object_ptr)
  end if
  this%cpp_object_ptr = c_null_ptr
end subroutine atlas_NextFunctionSpace__delete

subroutine atlas_NextFunctionSpace__reset(functionspace_out,functionspace_in)
  type(atlas_NextFunctionSpace), intent(inout) :: functionspace_out
  class(atlas_NextFunctionSpace), intent(in) :: functionspace_in
  if( .not. atlas_compare_equal(functionspace_out%cpp_object_ptr,functionspace_in%cpp_object_ptr) ) then
#ifndef FORTRAN_SUPPORTS_FINAL
    call atlas_NextFunctionSpace__finalize(functionspace_out)
#endif
    functionspace_out%cpp_object_ptr = functionspace_in%cpp_object_ptr
    if( c_associated(functionspace_out%cpp_object_ptr) ) call functionspace_out%attach()
  endif
end subroutine

function atlas_NextFunctionSpace__name(this) result(name)
  class(atlas_NextFunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__NextFunctionSpace__name(this%cpp_object_ptr)
  name = c_to_f_string_cptr(name_c_str)
end function


! -----------------------------------------------------------------------------
! FunctionSpace routines

function FunctionSpace__metadata(this) result(metadata)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  metadata%cpp_object_ptr = atlas__FunctionSpace__metadata(this%cpp_object_ptr)
end function FunctionSpace__metadata

subroutine FunctionSpace__create_field(this,name,nvars,kind)
  class(atlas_FunctionSpace), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: nvars
  integer, intent(in), optional :: kind
  if (present(kind)) then
    if (kind == ATLAS_KIND_REAL64) then
      call atlas__FunctionSpace__create_field_double(this%cpp_object_ptr,c_str(name),nvars)
    else if (kind == ATLAS_KIND_REAL32) then
      call atlas__FunctionSpace__create_field_float(this%cpp_object_ptr,c_str(name),nvars)
    else if (kind == ATLAS_KIND_INT32) then
      call atlas__FunctionSpace__create_field_int(this%cpp_object_ptr,c_str(name),nvars)
    else if (kind == ATLAS_KIND_INT64) then
      call atlas__FunctionSpace__create_field_long(this%cpp_object_ptr,c_str(name),nvars)
    else
      write(0,*) "Unsupported kind"
      write(0,*) 'call abort()'
    endif
  else if (wp == c_double) then
    call atlas__FunctionSpace__create_field_double(this%cpp_object_ptr,c_str(name),nvars)
  else if (wp == c_float) then
    call atlas__FunctionSpace__create_field_float(this%cpp_object_ptr,c_str(name),nvars)
  end if
end subroutine FunctionSpace__create_field


subroutine FunctionSpace__remove_field(this,name)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__FunctionSpace__remove_field(this%cpp_object_ptr,c_str(name))
end subroutine FunctionSpace__remove_field

function FunctionSpace__name(this) result(name)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__FunctionSpace__name(this%cpp_object_ptr)
  name = c_to_f_string_cptr(name_c_str)
end function FunctionSpace__name

function FunctionSpace__dof(this) result(dof)
  class(atlas_FunctionSpace), intent(in) :: this
  integer :: dof
  dof = atlas__FunctionSpace__dof(this%cpp_object_ptr)
end function FunctionSpace__dof

function FunctionSpace__glb_dof(this) result(glb_dof)
  class(atlas_FunctionSpace), intent(in) :: this
  integer :: glb_dof
  glb_dof = atlas__FunctionSpace__glb_dof(this%cpp_object_ptr)
end function FunctionSpace__glb_dof

function FunctionSpace__shape(this) result(shape)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, pointer :: shape(:)
  type(c_ptr) :: shape_c_ptr
  integer(c_int) :: field_rank
  call atlas__FunctionSpace__shapef(this%cpp_object_ptr, shape_c_ptr, field_rank)
  call C_F_POINTER ( shape_c_ptr , shape , (/field_rank/) )
end function FunctionSpace__shape


function FunctionSpace__field(this,name) result(field)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field = atlas_Field( atlas__FunctionSpace__field(this%cpp_object_ptr, c_str(name) ) )
  if( .not. C_associated(field%cpp_object_ptr) ) write(0,*) 'call abort()'
  call atlas_return(field)
end function FunctionSpace__field

function FunctionSpace__has_field(this,name) result(flag)
  class(atlas_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FunctionSpace__has_field(this%cpp_object_ptr, c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FunctionSpace__has_field

subroutine FunctionSpace__parallelise(this)
  class(atlas_FunctionSpace), intent(in) :: this
  call atlas__FunctionSpace__parallelise(this%cpp_object_ptr)
end subroutine FunctionSpace__parallelise

function FunctionSpace__get_halo_exchange(this) result(halo_exchange)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_HaloExchange) :: halo_exchange
  halo_exchange%cpp_object_ptr = atlas__FunctionSpace__halo_exchange( this%cpp_object_ptr )
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
  call atlas__FunctionSpace__halo_exchange_int( this%cpp_object_ptr, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r1
subroutine FunctionSpace__halo_exchange_int32_r2(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%cpp_object_ptr, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r2
subroutine FunctionSpace__halo_exchange_int32_r3(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%cpp_object_ptr, view, size(field_data) )
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
  call atlas__FunctionSpace__halo_exchange_float( this%cpp_object_ptr, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r1
subroutine FunctionSpace__halo_exchange_real32_r2(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%cpp_object_ptr, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r2
subroutine FunctionSpace__halo_exchange_real32_r3(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%cpp_object_ptr, view, size(field_data) )
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
  call atlas__FunctionSpace__halo_exchange_double( this%cpp_object_ptr, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r1
subroutine FunctionSpace__halo_exchange_real64_r2(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%cpp_object_ptr, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r2
subroutine FunctionSpace__halo_exchange_real64_r3(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%cpp_object_ptr, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r3
subroutine FunctionSpace__halo_exchange_real64_r4(this, field_data)
  class(atlas_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%cpp_object_ptr, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r4


function FunctionSpace__get_gather(this) result(gather)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_GatherScatter) :: gather
  gather%cpp_object_ptr = atlas__FunctionSpace__gather( this%cpp_object_ptr )
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
  call atlas__FunctionSpace__gather_float( this%cpp_object_ptr, field_data, size(field_data), &
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
  call atlas__FunctionSpace__gather_float( this%cpp_object_ptr, view, size(field_data), &
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
  call atlas__FunctionSpace__gather_float( this%cpp_object_ptr, view, size(field_data), &
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
  call atlas__FunctionSpace__gather_double( this%cpp_object_ptr, field_data, size(field_data), &
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
  call atlas__FunctionSpace__gather_double( this%cpp_object_ptr, view, size(field_data), &
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
  call atlas__FunctionSpace__gather_double( this%cpp_object_ptr, view, size(field_data), &
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
  call atlas__FunctionSpace__gather_int( this%cpp_object_ptr, field_data, size(field_data), &
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
  call atlas__FunctionSpace__gather_int( this%cpp_object_ptr, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r2

function FunctionSpace__get_checksum(this) result(checksum)
  class(atlas_FunctionSpace), intent(in) :: this
  type(atlas_Checksum) :: checksum
  checksum%cpp_object_ptr = atlas__FunctionSpace__checksum( this%cpp_object_ptr )
end function FunctionSpace__get_checksum
