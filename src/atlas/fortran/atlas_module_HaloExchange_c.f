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


subroutine HaloExchange__execute_int32_r1(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__HaloExchange__execute_int( this%private%object, field_data, 1 )
  else
    call atlas__HaloExchange__execute_int( this%private%object, field_data, nb_vars )
  end if
end subroutine HaloExchange__execute_int32_r1


subroutine HaloExchange__execute_int32_r2(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, intent(in) :: nb_vars
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_int( this%private%object, view, nb_vars )
end subroutine HaloExchange__execute_int32_r2


subroutine HaloExchange__execute_int32_r3(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, intent(in) :: nb_vars
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_int( this%private%object, view, nb_vars )
end subroutine HaloExchange__execute_int32_r3

subroutine HaloExchange__execute_real32_r1(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__HaloExchange__execute_float( this%private%object, field_data, 1 )
  else
    call atlas__HaloExchange__execute_float( this%private%object, field_data, nb_vars )
  end if
end subroutine HaloExchange__execute_real32_r1
subroutine HaloExchange__execute_real32_r2(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  integer, intent(in) :: nb_vars
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_float( this%private%object, view, nb_vars )
end subroutine HaloExchange__execute_real32_r2
subroutine HaloExchange__execute_real32_r3(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  integer, intent(in) :: nb_vars
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_float( this%private%object, view, nb_vars )
end subroutine HaloExchange__execute_real32_r3

subroutine HaloExchange__execute_real64_r1(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__HaloExchange__execute_double( this%private%object, field_data, 1 )
  else
    call atlas__HaloExchange__execute_double( this%private%object, field_data, nb_vars )
  end if
end subroutine HaloExchange__execute_real64_r1
subroutine HaloExchange__execute_real64_r2(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  integer, intent(in) :: nb_vars
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_double( this%private%object, view, nb_vars )
end subroutine HaloExchange__execute_real64_r2
subroutine HaloExchange__execute_real64_r3(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  integer, intent(in) :: nb_vars
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_double( this%private%object, view, nb_vars )
end subroutine HaloExchange__execute_real64_r3

! -----------------------------------------------------------------------------
