module atlas_C_interop
use, intrinsic :: iso_c_binding, only: c_int, c_long, c_float, c_double, c_ptr, &
  & c_loc, c_f_pointer, c_char, c_null_char, c_null_ptr
implicit none

private :: c_int, c_long, c_float, c_double, c_ptr, &
& c_loc, c_f_pointer, c_char, c_null_char, c_null_ptr

private


!========================================================================
! Public interface

!public atlas_free
!public atlas_compare_equal
!public get_c_arguments
!public c_to_f_string_str
!public c_to_f_string_cptr
!public c_str
!public c_str_no_trim
!public MAX_STR_LEN

!========================================================================

integer(c_int), parameter :: MAX_STR_LEN = 255
integer(c_int), target :: zero_length_array_int32(0)
integer(c_long),target :: zero_length_array_int64(0)
real(c_float),  target :: zero_length_array_real32(0)
real(c_double), target :: zero_length_array_real64(0)

interface
  subroutine atlas_free(ptr) bind(C)
    use, intrinsic :: iso_c_binding, only: c_ptr
    type(c_ptr), value :: ptr
  end subroutine atlas_free

  !int atlas__compare_cptr_equal( void* p1, void* p2 )
  function atlas__compare_cptr_equal(p1,p2) bind(c,name="atlas__compare_cptr_equal") result(equal)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    integer(c_int) :: equal
    type(c_ptr), value :: p1
    type(c_ptr), value :: p2
  end function

end interface

! =============================================================================
CONTAINS
! =============================================================================

function atlas_compare_equal(p1,p2) result(equal)
  use, intrinsic :: iso_c_binding, only: c_ptr
  logical :: equal
  type(c_ptr), intent(in) :: p1, p2
  if( atlas__compare_cptr_equal(p1,p2) == 1 ) then
    equal = .True.
  else
    equal = .False.
  endif
end function


subroutine get_c_arguments(argc,argv)
  integer(c_int), intent(out) :: argc
  type(c_ptr), intent(out) :: argv(:)
  character(kind=c_char,len=1), save, target :: args(255)
  character(kind=c_char,len=255), save, target :: cmd
  character(kind=c_char,len=255) :: arg
  integer :: iarg, arglen, pos, ich, argpos
  call get_command(cmd)
  do ich=1,len(cmd)
    if (cmd(ich:ich) == " ") then
      cmd(ich:ich) = c_null_char
      exit
    endif
  enddo
  argv(1) = c_loc(cmd(1:1))
  argc = command_argument_count()+1
  pos = 1
  do iarg=1,argc
    argpos = pos
    call get_command_argument(iarg, arg )
    arglen = len_trim(arg)
    do ich=1,arglen
      args(pos) = arg(ich:ich)
      pos = pos+1
    end do
    args(pos) = c_null_char;  pos = pos+1
    args(pos) = " ";          pos = pos+1
    argv(iarg+1) = c_loc(args(argpos))
  enddo
end subroutine


! -----------------------------------------------------------------------------
! Helper functions

function c_to_f_string_str(s) result(str)
  use, intrinsic :: iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  enddo
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  do i=1,nchars
    str(i:i) = s(i)
  enddo
end function c_to_f_string_str

function c_to_f_string_cptr(cptr) result(str)
  use, intrinsic :: iso_c_binding
  type(c_ptr), intent(in) :: cptr
  character(len=:), allocatable :: str
  character, dimension(:), pointer  :: s
  call C_F_POINTER ( cptr , s, (/MAX_STR_LEN/) )
  str = c_to_f_string_str(s)
end function c_to_f_string_cptr

function c_str(f_str)
  use, intrinsic :: iso_c_binding, only: c_char, c_null_char
  character(len=*), intent(in) :: f_str
  character(kind=c_char,len=len_trim(f_str)+1) :: c_str
  c_str = trim(f_str) // c_null_char
end function c_str

function c_str_no_trim(f_str)
  use, intrinsic :: iso_c_binding, only: c_char, c_null_char
  character(len=*), intent(in) :: f_str
  character(kind=c_char,len=len(f_str)+1) :: c_str_no_trim
  c_str_no_trim = f_str // c_null_char
end function c_str_no_trim

end module
