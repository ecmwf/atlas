module common_module
implicit none
save

private

integer, parameter, public :: jprb = selected_real_kind(13,300)
integer, parameter, public :: XX = 1
integer, parameter, public :: YY = 2

integer, parameter, public :: LOG_LEVEL_DEBUG   = 4
integer, parameter, public :: LOG_LEVEL_INFO    = 3
integer, parameter, public :: LOG_LEVEL_WARNING = 2
integer, parameter, public :: LOG_LEVEL_ERROR   = 1
integer, parameter, public :: LOG_UNIT = 0

public :: L2norm
public :: set_log_level
public :: log_error
public :: log_warning
public :: log_info
public :: log_debug

integer :: log_level = 3

public :: str

interface str
  module procedure str_integer
  module procedure str_real
end interface str


contains

function L2norm(array) result(norm)
  real(kind=jprb), intent(in) :: array(:)
  real(kind=jprb) :: norm
  integer :: i
  norm = 0
  do i=1,size(array)
    norm = norm + array(i)**2
  end do
  norm = sqrt(norm)/real(size(array))
end function L2norm

subroutine set_log_level(int)
  integer, intent(in) :: int
  log_level = int
end subroutine set_log_level

subroutine log_error(msg)
  character(len=*) :: msg
  if (log_level .ge. LOG_LEVEL_ERROR) then
    write(0,*) trim(msg)
  end if
end subroutine log_error

subroutine log_warning(msg)
  character(len=*) :: msg
  if (log_level .ge. LOG_LEVEL_WARNING) then
    write(0,*) trim(msg)
  end if
end subroutine log_warning

subroutine log_info(msg)
  character(len=*) :: msg
  if (log_level .ge. LOG_LEVEL_INFO) then
    write(0,*) trim(msg)
  end if
end subroutine log_info

subroutine log_debug(msg)
  character(len=*) :: msg
  if (log_level .ge. LOG_LEVEL_DEBUG) then
    write(0,*) trim(msg)
  end if
end subroutine log_debug

function str_integer(int) result(str)
  integer, intent(in) :: int
  character(len=10) :: str
  write(str,*) int
end function str_integer

function str_real(re) result(str)
  real(kind=jprb), intent(in) :: re
  character(len=20) :: str
  write(str,*) re
end function str_real

end module common_module