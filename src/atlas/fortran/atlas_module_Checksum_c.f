! (C) Copyright 2013-2014 ECMWF.

! ------------------------------------------------------------------------------
! Checksum routines

function new_Checksum() result(Checksum)
  type(Checksum_type) :: Checksum
  Checksum%private%object = atlas__Checksum__new()
end function new_Checksum

subroutine Checksum__delete(this)
  type(Checksum_type), intent(inout) :: this
  if ( c_associated(this%private%object) ) then
    call atlas__Checksum__delete(this%private%object)
  end if
  this%private%object = C_NULL_ptr
end subroutine Checksum__delete

subroutine Checksum__setup(this, part, remote_idx, glb_idx, opt_max_glb_idx)
  class(Checksum_type), intent(in) :: this
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
  call atlas__Checksum__setup( this%private%object, part, remote_idx, 1, &
    &                          glb_idx, max_glb_idx, size(part) )
end subroutine Checksum__setup

function Checksum__execute_int32_r1(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer :: lstrides(1), lextents(1), lrank=1
  lstrides = (/ stride(loc_field_data,2) /)
  lextents = (/ 1                        /)
  call atlas__Checksum__execute_strided_int( this%private%object, &
    &  loc_field_data, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_int32_r1

function Checksum__execute_int32_r2(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer, pointer :: lview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  lstrides = (/ stride(loc_field_data,2), stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_int( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_int32_r2


function Checksum__execute_int32_r3(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer, pointer :: lview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  lstrides = (/ stride(loc_field_data,3), stride(loc_field_data,2) , stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size(loc_field_data,1) /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_int( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_int32_r3

function Checksum__execute_real32_r1(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  real(c_float), intent(in)   :: loc_field_data(:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer :: lstrides(1), lextents(1), lrank=1
  lstrides = (/ stride(loc_field_data,1) /)
  lextents = (/ 1                        /)
  call atlas__Checksum__execute_strided_float( this%private%object, &
    &  loc_field_data, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_real32_r1
function Checksum__execute_real32_r2(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_float), pointer :: lview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  lstrides = (/ stride(loc_field_data,2), stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_float( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_real32_r2
function Checksum__execute_real32_r3(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_float), pointer :: lview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  lstrides = (/ stride(loc_field_data,3), stride(loc_field_data,2) , stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_float( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_real32_r3

function Checksum__execute_real64_r1(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  real(c_double), intent(in)   :: loc_field_data(:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer :: lstrides(1), lextents(1), lrank=1
  real(c_double), pointer :: lview(:)
  lstrides = (/ stride(loc_field_data,1) /)
  lextents = (/ 1                        /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_double( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_real64_r1
function Checksum__execute_real64_r2(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_double), pointer :: lview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  lstrides = (/ stride(loc_field_data,2), stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_double( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_real64_r2
function Checksum__execute_real64_r3(this, loc_field_data) result(checksum)
  class(Checksum_type), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_double), pointer :: lview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  lstrides = (/ stride(loc_field_data,3), stride(loc_field_data,2) , stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size  (loc_field_data,1) /)
  lview => view1d(loc_field_data)
  call atlas__Checksum__execute_strided_double( this%private%object, &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_to_f_string_str(checksum_c_str)
end function Checksum__execute_real64_r3

! -----------------------------------------------------------------------------
