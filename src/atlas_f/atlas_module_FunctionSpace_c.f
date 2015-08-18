! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
!

function Nodes__metadata(this) result(metadata)
  use atlas_nodes_c_binding
  type(atlas_Metadata) :: metadata
  class(atlas_Nodes), intent(in) :: this
  metadata%cpp_object_ptr = atlas__Nodes__metadata(this%cpp_object_ptr)
end function

function Nodes__size(this) result(size)
  use atlas_nodes_c_binding
  integer :: size
  class(atlas_Nodes), intent(in) :: this
  size = atlas__Nodes__size(this%cpp_object_ptr)
end function

subroutine Nodes__resize(this,size)
  use atlas_nodes_c_binding
  class(atlas_Nodes), intent(in) :: this
  integer, intent(in) :: size
  call atlas__Nodes__resize(this%cpp_object_ptr,size)
end subroutine

function Nodes__add(this,field) result(added_field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: added_field
  class(atlas_Nodes), intent(inout) :: this
  type(atlas_Field), intent(in) :: field
  added_field%cpp_object_ptr = atlas__Nodes__add(this%cpp_object_ptr,field%cpp_object_ptr)
end function

subroutine Nodes__remove_field(this,name)
  use atlas_nodes_c_binding
  class(atlas_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__Nodes__remove_field(this%cpp_object_ptr,c_str(name))
end subroutine

function Nodes__has_field(this,name) result(has_field)
  use atlas_nodes_c_binding
  logical :: has_field
  class(atlas_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  integer :: has_field_int
  has_field_int = atlas__Nodes__has_field(this%cpp_object_ptr,c_str(name))
  has_field = .False.
  if( has_field_int == 1 ) has_field = .True.
end function

function Nodes__nb_fields(this) result(nb_fields)
  use atlas_nodes_c_binding
  integer :: nb_fields
  class(atlas_Nodes), intent(in) :: this
  nb_fields = atlas__Nodes__nb_fields(this%cpp_object_ptr)
end function

function Nodes__field_by_name(this,name) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  field%cpp_object_ptr = atlas__Nodes__field_by_name(this%cpp_object_ptr,c_str(name))
end function

function Nodes__field_by_idx(this,idx) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_Nodes), intent(in) :: this
  integer, intent(in) :: idx
  field%cpp_object_ptr = atlas__Nodes__field_by_idx(this%cpp_object_ptr,idx)
end function

function Nodes__str(this) result(str)
  use atlas_nodes_c_binding
  character(len=:), allocatable :: str
  class(atlas_Nodes), intent(in) :: this
  type(c_ptr) :: str_cptr
  integer(c_int) :: str_size
  call atlas__Nodes__str(this%cpp_object_ptr,str_cptr,str_size)
  allocate(character(len=str_size) :: str )
  str = c_to_f_string_cptr(str_cptr)
  call atlas_free(str_cptr)
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
  field%cpp_object_ptr = atlas__FunctionSpace__field(this%cpp_object_ptr, c_str(name) )
  if( .not. C_associated(field%cpp_object_ptr) ) write(0,*) 'call abort()'
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
