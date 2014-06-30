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

subroutine Gather__execute_int32_r1_r1(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:)
  integer, intent(out) :: glb_field_data(:)
  integer :: lstrides(1), lextents(1), lrank=1
  integer :: gstrides(1), gextents(1), grank=1
  integer(c_int), pointer :: lview(:), gview(:)
  lstrides = (/ stride(loc_field_data,2) /)
  lextents = (/ 1                        /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,2) /)
  gextents = (/ 1                        /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_int( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_int32_r1_r1


subroutine Gather__execute_int32_r2_r2(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:)
  integer, intent(out) :: glb_field_data(:,:)
  integer, pointer :: lview(:), gview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  integer :: gstrides(2), gextents(2), grank=2
  lstrides = (/ stride(loc_field_data,2), stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,2), stride(glb_field_data,1) /)
  gextents = (/ 1,                        size  (glb_field_data,1) /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_int( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_int32_r2_r2


subroutine Gather__execute_int32_r3_r3(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:,:)
  integer, intent(out) :: glb_field_data(:,:,:)
  integer, pointer :: lview(:), gview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  integer :: gstrides(3), gextents(3), grank=3
  lstrides = (/ stride(loc_field_data,3), stride(loc_field_data,2) , stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size(loc_field_data,1) /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,3), stride(glb_field_data,2) , stride(glb_field_data,1) /)
  gextents = (/ 1,                        size  (glb_field_data,2) , size  (glb_field_data,1) /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_int( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_int32_r3_r3

subroutine Gather__execute_real32_r1_r1(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:)
  real(c_float), intent(out) :: glb_field_data(:)
  integer :: lstrides(1), lextents(1), lrank=1
  integer :: gstrides(1), gextents(1), grank=1
  real(c_float), pointer :: lview(:), gview(:)
  lstrides = (/ stride(loc_field_data,2) /)
  lextents = (/ 1                        /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,2) /)
  gextents = (/ 1                        /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_float( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_real32_r1_r1
subroutine Gather__execute_real32_r2_r2(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:)
  real(c_float), intent(out) :: glb_field_data(:,:)
  real(c_float), pointer :: lview(:), gview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  integer :: gstrides(2), gextents(2), grank=2
  lstrides = (/ stride(loc_field_data,2), stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,2), stride(glb_field_data,1) /)
  gextents = (/ 1,                        size  (glb_field_data,1) /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_float( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_real32_r2_r2
subroutine Gather__execute_real32_r3_r3(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:,:)
  real(c_float), intent(out) :: glb_field_data(:,:,:)
  real(c_float), pointer :: lview(:), gview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  integer :: gstrides(3), gextents(3), grank=3
  lstrides = (/ stride(loc_field_data,3), stride(loc_field_data,2) , stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,3), stride(glb_field_data,2) , stride(glb_field_data,1) /)
  gextents = (/ 1,                        size  (glb_field_data,2) , size  (glb_field_data,1) /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_float( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_real32_r3_r3

subroutine Gather__execute_real64_r1_r1(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  real(c_double), intent(in)   :: loc_field_data(:)
  real(c_double), intent(out)  :: glb_field_data(:)
  integer :: lstrides(1), lextents(1), lrank=1
  integer :: gstrides(1), gextents(1), grank=1
  real(c_double), pointer :: lview(:), gview(:)
  lstrides = (/ stride(loc_field_data,1) /)
  lextents = (/ 1                        /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,1) /)
  gextents = (/ 1                        /)
  gview => view1d(glb_field_data)
!  write(0,*) MPL_rank(),"lstrides",lstrides
!  write(0,*) MPL_rank(),"lextents",lextents
!  write(0,*) MPL_rank(),"gstrides",gstrides
!  write(0,*) MPL_rank(),"gextents",gextents
!  write(0,*) MPL_rank(),"localsize",lstrides(1)*lextents(1)*size(loc_field_data)
!  write(0,*) "address, size = ",loc(loc_field_data(1)),size(loc_field_data), loc(lview(1))

  call atlas__Gather__execute_strided_double( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_real64_r1_r1
subroutine Gather__execute_real64_r2_r2(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:)
  real(c_double), intent(out) :: glb_field_data(:,:)
  real(c_double), pointer :: lview(:), gview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  integer :: gstrides(2), gextents(2), grank=2
  lstrides = (/ stride(loc_field_data,2), stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,2), stride(glb_field_data,1) /)
  gextents = (/ 1,                        size  (glb_field_data,1) /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_double( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_real64_r2_r2
subroutine Gather__execute_real64_r3_r3(this, loc_field_data, glb_field_data)
  class(Gather_type), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:,:)
  real(c_double), intent(out) :: glb_field_data(:,:,:)
  real(c_double), pointer :: lview(:), gview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  integer :: gstrides(3), gextents(3), grank=3
  lstrides = (/ stride(loc_field_data,3), stride(loc_field_data,2) , stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  gstrides = (/ stride(glb_field_data,3), stride(glb_field_data,2) , stride(glb_field_data,1) /)
  gextents = (/ 1,                        size  (glb_field_data,2) , size  (glb_field_data,1) /)
  gview => view1d(glb_field_data)
  if( size(gview) == 0 ) then
    allocate(gview(0))
  endif
  call atlas__Gather__execute_strided_double( this%private%object, &
    &  lview, lstrides, lextents, lrank, &
    &  gview, gstrides, gextents, grank )
end subroutine Gather__execute_real64_r3_r3

! -----------------------------------------------------------------------------
