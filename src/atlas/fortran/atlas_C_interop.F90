module atlas_C_interop
use, intrinsic :: iso_c_binding
implicit none

integer, private, parameter :: MAX_STR_LEN = 255
integer, private, parameter :: FIELD_NB_VARS = -1

INTERFACE view1d
  module procedure view1d_int32_rank2
  module procedure view1d_int32_rank3
  module procedure view1d_real32_rank2
  module procedure view1d_real32_rank3
  module procedure view1d_real64_rank2
  module procedure view1d_real64_rank3
  module procedure view1d_real64_rank4
end interface view1d

INTERFACE stride
  module procedure stride_int32_r1
  module procedure stride_int32_r2
  module procedure stride_int32_r3
  module procedure stride_real32_r1
  module procedure stride_real32_r2
  module procedure stride_real32_r3
  module procedure stride_real64_r1
  module procedure stride_real64_r2
  module procedure stride_real64_r3
end interface stride

! =============================================================================
CONTAINS
! =============================================================================




! -----------------------------------------------------------------------------
! Helper functions

function c_to_f_string_str(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  character(len=:), allocatable :: mold
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  allocate(character(len=nchars) :: mold)
  str = transfer(s(1:nchars), mold)
end function c_to_f_string_str

function c_to_f_string_cptr(cptr) result(str)
  use iso_c_binding
  type(c_ptr), intent(in) :: cptr
  character(len=:), allocatable :: str
  character, dimension(:), pointer  :: s
  call C_F_POINTER ( cptr , s, (/MAX_STR_LEN/) )
  str = c_to_f_string_str(s)
end function c_to_f_string_cptr

function c_str(f_str)
  use, intrinsic :: iso_c_binding
  character(len=*), intent(in) :: f_str
  character(len=len_trim(f_str)+1) :: c_str
  c_str = trim(f_str) // c_null_char
end function c_str

! ------------------------------------------------------------
! view interface

function c_loc_int32(x)
  use iso_c_binding
  integer, target :: x
  type(c_ptr) :: c_loc_int32
  c_loc_int32 = C_LOC(x)
end function

function c_loc_real32(x)
  use iso_c_binding
  real(c_float), target :: x
  type(c_ptr) :: c_loc_real32
  c_loc_real32 = C_LOC(x)
end function

function c_loc_real64(x)
  use iso_c_binding
  real(c_double), target :: x
  type(c_ptr) :: c_loc_real64
  c_loc_real64 = C_LOC(x)
end function

function view1d_int32_rank2(array) result( view )
  integer, intent(in), target :: array(:,:)
  type(c_ptr) :: array_c_ptr
  integer, pointer :: view(:)
  array_c_ptr = c_loc_int32(array(1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_int32_rank2

function view1d_int32_rank3(array) result( view )
  integer, intent(in), target :: array(:,:,:)
  type(c_ptr) :: array_c_ptr
  integer, pointer :: view(:)
  array_c_ptr = c_loc_int32(array(1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_int32_rank3

function view1d_real32_rank2(array) result( view )
  real(c_float), intent(in), target :: array(:,:)
  type(c_ptr) :: array_c_ptr
  real(c_float), pointer :: view(:)
#ifndef  __GFORTRAN__
  if( .not. is_contiguous(array) ) then
    write(0,*) "ERROR: array is not contiguous in view1d"
    write(0,*) 'call abort()'
  end if
#endif
  array_c_ptr = c_loc_real32(array(1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real32_rank2

function view1d_real32_rank3(array) result( view )
  real(c_float), intent(in), target :: array(:,:,:)
  type(c_ptr) :: array_c_ptr
  real(c_float), pointer :: view(:)
  array_c_ptr = c_loc_real32(array(1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real32_rank3

function view1d_real64_rank2(array) result( view )
  real(c_double), intent(in), target :: array(:,:)
  type(c_ptr) :: array_c_ptr
  real(c_double), pointer :: view(:)
  array_c_ptr = c_loc_real64(array(1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real64_rank2

function view1d_real64_rank3(array) result( view )
  real(c_double), intent(in), target :: array(:,:,:)
  type(c_ptr) :: array_c_ptr
  real(c_double), pointer :: view(:)
  array_c_ptr = c_loc_real64(array(1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real64_rank3

function view1d_real64_rank4(array) result( view )
  real(c_double), intent(in), target :: array(:,:,:,:)
  type(c_ptr) :: array_c_ptr
  real(c_double), pointer :: view(:)
  array_c_ptr = c_loc_real64(array(1,1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real64_rank4

! ------------------------------------------------------------
! stride interface

function stride_int32_r1(arr,dim) result( stride )
  integer :: arr(:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2))-loc(arr(1)))/4
end function stride_int32_r1

function stride_int32_r2(arr,dim) result( stride )
  integer :: arr(:,:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2,1))-loc(arr(1,1)))/4
  if (dim == 2) stride = (loc(arr(1,2))-loc(arr(1,1)))/4
end function stride_int32_r2

function stride_int32_r3(arr,dim) result( stride )
  integer :: arr(:,:,:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2,1,1))-loc(arr(1,1,1)))/4
  if (dim == 2) stride = (loc(arr(1,2,1))-loc(arr(1,1,1)))/4
  if (dim == 3) stride = (loc(arr(1,1,2))-loc(arr(1,1,1)))/4
end function stride_int32_r3

function stride_real32_r1(arr,dim) result( stride )
  real(c_float) :: arr(:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2))-loc(arr(1)))/4
end function stride_real32_r1

function stride_real32_r2(arr,dim) result( stride )
  real(c_float) :: arr(:,:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2,1))-loc(arr(1,1)))/4
  if (dim == 2) stride = (loc(arr(1,2))-loc(arr(1,1)))/4
end function stride_real32_r2

function stride_real32_r3(arr,dim) result( stride )
  real(c_float) :: arr(:,:,:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2,1,1))-loc(arr(1,1,1)))/4
  if (dim == 2) stride = (loc(arr(1,2,1))-loc(arr(1,1,1)))/4
  if (dim == 3) stride = (loc(arr(1,1,2))-loc(arr(1,1,1)))/4
end function stride_real32_r3

function stride_real64_r1(arr,dim) result( stride )
  real(c_double) :: arr(:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2))-loc(arr(1)))/8
end function stride_real64_r1

function stride_real64_r2(arr,dim) result( stride )
  real(c_double) :: arr(:,:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2,1))-loc(arr(1,1)))/8
  if (dim == 2) stride = (loc(arr(1,2))-loc(arr(1,1)))/8
end function stride_real64_r2

function stride_real64_r3(arr,dim) result( stride )
  real(c_double) :: arr(:,:,:)
  integer :: dim
  integer :: stride
  if (dim == 1) stride = (loc(arr(2,1,1))-loc(arr(1,1,1)))/8
  if (dim == 2) stride = (loc(arr(1,2,1))-loc(arr(1,1,1)))/8
  if (dim == 3) stride = (loc(arr(1,1,2))-loc(arr(1,1,1)))/8
end function stride_real64_r3

end module
