! (C) Copyright 2013-2014 ECMWF.

! ------------------------------------------------------------------------------
! Gather routines

function new_Gather() result(gather)
  type(Gather_type) :: gather
  gather%private%object = atlas__Gather__new()
end function new_Gather

subroutine Gather__delete(this)
  type(Gather_type), intent(inout) :: this
  if ( c_associated(this%private%object) ) then
    call atlas__Gather__delete(this%private%object)
  end if
  this%private%object = C_NULL_ptr
end subroutine Gather__delete

subroutine Gather__setup(this, part, remote_idx, glb_idx, opt_max_glb_idx)
  class(Gather_type), intent(in) :: this
  integer, intent(in) :: part(:)
  integer, intent(in) :: remote_idx(:)
  integer, intent(in) :: glb_idx(:)
  integer, optional, intent(in) :: opt_max_glb_idx
  integer :: max_glb_idx
  if (.not. present(opt_max_glb_idx) ) then
    max_glb_idx = huge(max_glb_idx)
  else
    max_glb_idx = opt_max_glb_idx
  endif
  call atlas__Gather__setup( this%private%object, part, remote_idx, 1, &
    &                        glb_idx, max_glb_idx, size(part) )
end subroutine Gather__setup

function Gather__glb_dof(this) result(glb_dof)
  class(Gather_type), intent(in) :: this
  integer :: glb_dof
  glb_dof = atlas__Gather__glb_dof(this%private%object)
end function Gather__glb_dof

subroutine Gather__execute_int32_r1(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:)
  integer, intent(out) :: glb_field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__Gather__execute_int( this%private%object, loc_field_data, glb_field_data, 1 )
  else
    call atlas__Gather__execute_int( this%private%object, loc_field_data, glb_field_data, nb_vars )
  end if
end subroutine Gather__execute_int32_r1


subroutine Gather__execute_int32_r2(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:)
  integer, intent(out) :: glb_field_data(:,:)
  integer, intent(in)  :: nb_vars
  integer, pointer :: lview(:), gview(:)
  lview => view1d(loc_field_data)
  gview => view1d(glb_field_data)
  call atlas__Gather__execute_int( this%private%object, lview, gview, nb_vars )
end subroutine Gather__execute_int32_r2


subroutine Gather__execute_int32_r3(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:,:)
  integer, intent(out) :: glb_field_data(:,:,:)
  integer, intent(in)  :: nb_vars
  integer, pointer :: lview(:), gview(:)
  lview => view1d(loc_field_data)
  gview => view1d(glb_field_data)
  call atlas__Gather__execute_int( this%private%object, lview, gview, nb_vars )
end subroutine Gather__execute_int32_r3

subroutine Gather__execute_real32_r1(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:)
  real(c_float), intent(out) :: glb_field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__Gather__execute_float( this%private%object, loc_field_data, glb_field_data, 1 )
  else
    call atlas__Gather__execute_float( this%private%object, loc_field_data, glb_field_data, nb_vars )
  end if
end subroutine Gather__execute_real32_r1
subroutine Gather__execute_real32_r2(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:)
  real(c_float), intent(out) :: glb_field_data(:,:)
  integer, intent(in) :: nb_vars
  real(c_float), pointer :: lview(:), gview(:)
  lview => view1d(loc_field_data)
  gview => view1d(glb_field_data)
  call atlas__Gather__execute_float( this%private%object, lview, gview, nb_vars )
end subroutine Gather__execute_real32_r2
subroutine Gather__execute_real32_r3(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:,:)
  real(c_float), intent(out) :: glb_field_data(:,:,:)
  integer, intent(in) :: nb_vars
  real(c_float), pointer :: lview(:), gview(:)
  lview => view1d(loc_field_data)
  gview => view1d(glb_field_data)
  call atlas__Gather__execute_float( this%private%object, lview, gview, nb_vars )
end subroutine Gather__execute_real32_r3

subroutine Gather__execute_real64_r1(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  real(c_double), intent(in)   :: loc_field_data(:)
  real(c_double), intent(out)  :: glb_field_data(:)
  integer, optional, intent(in)  :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__Gather__execute_double( this%private%object, loc_field_data, glb_field_data, 1 )
  else
    call atlas__Gather__execute_double( this%private%object, loc_field_data, glb_field_data, nb_vars )
  end if
end subroutine Gather__execute_real64_r1
subroutine Gather__execute_real64_r2(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:)
  real(c_double), intent(out) :: glb_field_data(:,:)
  integer, intent(in) :: nb_vars
  real(c_double), pointer :: lview(:), gview(:)
  lview => view1d(loc_field_data)
  gview => view1d(glb_field_data)
  call atlas__Gather__execute_double( this%private%object, lview, gview, nb_vars )
end subroutine Gather__execute_real64_r2
subroutine Gather__execute_real64_r3(this, loc_field_data, glb_field_data, nb_vars)
  class(Gather_type), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:,:)
  real(c_double), intent(out) :: glb_field_data(:,:,:)
  integer, intent(in) :: nb_vars
  real(c_double), pointer :: lview(:), gview(:)
  lview => view1d(loc_field_data)
  gview => view1d(glb_field_data)
  call atlas__Gather__execute_double( this%private%object, lview, gview, nb_vars )
end subroutine Gather__execute_real64_r3

! -----------------------------------------------------------------------------
