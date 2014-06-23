! (C) Copyright 2013-2014 ECMWF.

! ------------------------------------------------------------------------------
! HaloExchange routines

function new_HaloExchange() result(halo_exchange)
  type(HaloExchange_type) :: halo_exchange
  halo_exchange%private%object = atlas__HaloExchange__new()
end function new_HaloExchange

subroutine HaloExchange__delete(this)
  type(HaloExchange_type), intent(inout) :: this
  if ( c_associated(this%private%object) ) then
    call atlas__HaloExchange__delete(this%private%object)
  end if
  this%private%object = C_NULL_ptr
end subroutine HaloExchange__delete

subroutine HaloExchange__setup(this, part, remote_idx)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(in) :: part(:)
  integer, intent(in) :: remote_idx(:)
  call atlas__HaloExchange__setup( this%private%object, part, remote_idx, 1, size(part) )
end subroutine HaloExchange__setup


subroutine HaloExchange__execute_int32_r1(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_int( this%private%object, field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_int32_r1

subroutine HaloExchange__execute_int32_r2(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_int( this%private%object, view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_int32_r2

subroutine HaloExchange__execute_int32_r3(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  integer :: strides(3), extents(3)
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_int( this%private%object, view, &
    & strides, extents, 3 )
end subroutine HaloExchange__execute_int32_r3

subroutine HaloExchange__execute_real32_r1(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_float( this%private%object, field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_real32_r1
subroutine HaloExchange__execute_real32_r2(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_float( this%private%object, view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_real32_r2
subroutine HaloExchange__execute_real32_r3(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  integer :: strides(3), extents(3), rank=3
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_float( this%private%object, view, &
    & strides, extents, rank )
end subroutine HaloExchange__execute_real32_r3

subroutine HaloExchange__execute_real64_r1(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_double( this%private%object, field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_real64_r1
subroutine HaloExchange__execute_real64_r2(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_double( this%private%object, view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_real64_r2
subroutine HaloExchange__execute_real64_r3(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  integer :: strides(3), extents(3), rank=3
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_double( this%private%object, view, &
    & strides, extents, rank )
end subroutine HaloExchange__execute_real64_r3
subroutine HaloExchange__execute_real64_r4(this, field_data)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  integer :: strides(4), extents(4), rank=4
  view => view1d(field_data)
  strides = (/ stride(field_data,4), stride(field_data,3), stride(field_data,2), stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,3), ubound(field_data,2), ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_double( this%private%object, view, &
    & strides, extents, rank )
end subroutine HaloExchange__execute_real64_r4

! -----------------------------------------------------------------------------
