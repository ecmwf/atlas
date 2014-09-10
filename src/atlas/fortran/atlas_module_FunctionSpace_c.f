! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! FunctionSpace routines

function new_FunctionSpace(name,shape_func,nb_nodes) result(function_space)
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_nodes
  type(FunctionSpace_type) :: function_space
  integer :: extents(2)
  extents = (/nb_nodes,FIELD_NB_VARS/)
  function_space%private%object = atlas__FunctionSpace__new(c_str(name),c_str(shape_func), &
    & extents, 2 )
end function new_FunctionSpace

function new_PrismaticFunctionSpace(name,shape_func,nb_levels,nb_nodes) result(function_space)
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_levels
  integer, intent(in) :: nb_nodes
  type(FunctionSpace_type) :: function_space
  integer :: extents(3)
  extents = (/nb_nodes,nb_levels,FIELD_NB_VARS/)
  function_space%private%object = atlas__FunctionSpace__new(c_str(name),c_str(shape_func), &
    & extents, 3 )
end function new_PrismaticFunctionSpace

subroutine FunctionSpace__delete(this)
  type(FunctionSpace_type), intent(inout) :: this
  if ( c_associated(this%private%object) ) then
    call atlas__FunctionSpace__delete(this%private%object)
  end if
  this%private%object = C_NULL_ptr
end subroutine FunctionSpace__delete

subroutine FunctionSpace__create_field(this,name,nvars,kind)
  class(FunctionSpace_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: nvars
  integer, intent(in), optional :: kind
  if (present(kind)) then
    if (kind == KIND_REAL64) then
      call atlas__FunctionSpace__create_field_double(this%private%object,c_str(name),nvars)
    else if (kind == KIND_REAL32) then
      call atlas__FunctionSpace__create_field_float(this%private%object,c_str(name),nvars)
    else if (kind == KIND_INT32) then
      call atlas__FunctionSpace__create_field_int(this%private%object,c_str(name),nvars)
    else
      write(0,*) "Unsupported kind"
      write(0,*) 'call abort()'
    endif
  else if (wp == c_double) then
    call atlas__FunctionSpace__create_field_double(this%private%object,c_str(name),nvars)
  else if (wp == c_float) then
    call atlas__FunctionSpace__create_field_float(this%private%object,c_str(name),nvars)
  end if
end subroutine FunctionSpace__create_field


subroutine FunctionSpace__remove_field(this,name)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__FunctionSpace__remove_field(this%private%object,c_str(name))
end subroutine FunctionSpace__remove_field

function FunctionSpace__name(this) result(name)
  class(FunctionSpace_type), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__FunctionSpace__name(this%private%object)
  name = c_to_f_string_cptr(name_c_str)
end function FunctionSpace__name

function FunctionSpace__dof(this) result(dof)
  class(FunctionSpace_type), intent(in) :: this
  integer :: dof
  dof = atlas__FunctionSpace__dof(this%private%object)
end function FunctionSpace__dof

function FunctionSpace__glb_dof(this) result(glb_dof)
  class(FunctionSpace_type), intent(in) :: this
  integer :: glb_dof
  glb_dof = atlas__FunctionSpace__glb_dof(this%private%object)
end function FunctionSpace__glb_dof

function FunctionSpace__bounds(this) result(bounds)
  class(FunctionSpace_type), intent(in) :: this
  integer, pointer :: bounds(:)
  type(c_ptr) :: bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__FunctionSpace__shapef(this%private%object, bounds_c_ptr, field_rank)
  call C_F_POINTER ( bounds_c_ptr , bounds , (/field_rank/) )
end function FunctionSpace__bounds


function FunctionSpace__field(this,name) result(field)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(Field_type) :: field
  field%private%object = atlas__FunctionSpace__field(this%private%object, c_str(name) )
  if( .not. C_associated(field%private%object) ) write(0,*) 'call abort()'
end function FunctionSpace__field

function FunctionSpace__has_field(this,name) result(flag)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FunctionSpace__has_field(this%private%object, c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FunctionSpace__has_field

subroutine FunctionSpace__parallelise(this)
  class(FunctionSpace_type), intent(in) :: this
  call atlas__FunctionSpace__parallelise(this%private%object)
end subroutine FunctionSpace__parallelise

function FunctionSpace__get_halo_exchange(this) result(halo_exchange)
  class(FunctionSpace_type), intent(in) :: this
  type(HaloExchange_type) :: halo_exchange
  halo_exchange%private%object = atlas__FunctionSpace__halo_exchange( this%private%object )
end function FunctionSpace__get_halo_exchange

subroutine FunctionSpace__halo_exchange_int32_r1(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_int( this%private%object, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r1
subroutine FunctionSpace__halo_exchange_int32_r2(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r2
subroutine FunctionSpace__halo_exchange_int32_r3(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r3

subroutine FunctionSpace__halo_exchange_real32_r1(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_float( this%private%object, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r1
subroutine FunctionSpace__halo_exchange_real32_r2(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r2
subroutine FunctionSpace__halo_exchange_real32_r3(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r3

subroutine FunctionSpace__halo_exchange_real64_r1(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_double( this%private%object, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r1
subroutine FunctionSpace__halo_exchange_real64_r2(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r2
subroutine FunctionSpace__halo_exchange_real64_r3(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r3
subroutine FunctionSpace__halo_exchange_real64_r4(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%private%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r4


function FunctionSpace__get_gather(this) result(gather)
  class(FunctionSpace_type), intent(in) :: this
  type(GatherScatter_type) :: gather
  gather%private%object = atlas__FunctionSpace__gather( this%private%object )
end function FunctionSpace__get_gather

subroutine FunctionSpace__gather_real32_r1(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(in) :: field_data(:)
  real(c_float), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_float( this%private%object, field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r1


subroutine FunctionSpace__gather_real32_r2(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_float( this%private%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r2


subroutine FunctionSpace__gather_real32_r3(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_float( this%private%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r3

subroutine FunctionSpace__gather_real64_r1(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(in) :: field_data(:)
  real(c_double), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_double( this%private%object, field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r1


subroutine FunctionSpace__gather_real64_r2(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_double( this%private%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r2


subroutine FunctionSpace__gather_real64_r3(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_double( this%private%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r3

subroutine FunctionSpace__gather_int32_r1(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(in) :: field_data(:)
  integer, intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_int( this%private%object, field_data, size(field_data), &
                                       & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r1

subroutine FunctionSpace__gather_int32_r2(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(in) :: field_data(:,:)
  integer, intent(inout) :: glbfield_data(:,:)
  integer, pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_int( this%private%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r2

function FunctionSpace__get_checksum(this) result(checksum)
  class(FunctionSpace_type), intent(in) :: this
  type(Checksum_type) :: checksum
  checksum%private%object = atlas__FunctionSpace__checksum( this%private%object )
end function FunctionSpace__get_checksum
