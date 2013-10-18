module common_module

use parallel_module, only: myproc

implicit none
save

private

integer, parameter, public :: jprb = selected_real_kind(13,300)
integer, parameter, public :: jpim = selected_int_kind(9)

integer, parameter, public :: XX = 1
integer, parameter, public :: YY = 2

integer, parameter, public :: LOG_LEVEL_DEBUG   = 4
integer, parameter, public :: LOG_LEVEL_INFO    = 3
integer, parameter, public :: LOG_LEVEL_WARNING = 2
integer, parameter, public :: LOG_LEVEL_ERROR   = 1
integer, parameter, public :: LOG_UNIT = 0

public :: L2norm
public :: set_log_level
public :: set_log_proc
public :: log_error
public :: log_warning
public :: log_info
public :: log_debug
public :: Timer_type

integer :: log_level = 3
integer :: log_proc = -1 ! means everyone
character(len=1024), public :: log_str
public :: str

interface str
  module procedure str_integer
  module procedure str_real
end interface str

type :: Timer_type
private
  integer :: clck_counts_start, clck_counts_stop, clck_rate
  integer :: counted = 0
  logical :: paused = .True.
contains
  procedure, public :: start   => Timer_start
  procedure, public :: pause   => Timer_pause
  procedure, public :: resume  => Timer_resume
  procedure, public :: elapsed => Timer_elapsed
end type Timer_type


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

function Timer_elapsed(self) result(time)
  class(Timer_type), intent(inout) :: self
  real(kind=jprb) :: time
  if (.not. self%paused) then
    call system_clock ( self%clck_counts_stop, self%clck_rate )
    time = (self%counted + self%clck_counts_stop - self%clck_counts_start)/real(self%clck_rate)
  else if (self%counted .ge. 0) then
    time = self%counted/real(self%clck_rate)
  else
    time = 0.
  end if
end function Timer_elapsed

subroutine Timer_start(self)
  class(Timer_type), intent(inout) :: self
  call system_clock ( self%clck_counts_start, self%clck_rate )
  self%paused = .False.
  self%counted = 0
end subroutine Timer_start

subroutine Timer_pause(self)
  class(Timer_type), intent(inout) :: self
  call system_clock ( self%clck_counts_stop, self%clck_rate )
  self%counted = self%counted + self%clck_counts_stop - self%clck_counts_start
  self%paused = .True.
end subroutine Timer_pause

subroutine Timer_resume(self)
  class(Timer_type), intent(inout) :: self
  call system_clock ( self%clck_counts_start, self%clck_rate )
  self%paused = .False.
end subroutine Timer_resume

subroutine set_log_level(int)
  integer, intent(in) :: int
  log_level = int
end subroutine set_log_level

subroutine set_log_proc(int)
  integer, intent(in) :: int
  log_proc = int
end subroutine set_log_proc

subroutine log_error(msg)
  character(len=*), optional :: msg
  if (log_level .ge. LOG_LEVEL_ERROR) then
    if( log_proc < 0 .or. log_proc .eq. myproc ) then
      if (present(msg)) then
        write(0,'(A)') trim(msg)
      else 
        write(0,'(A)') trim(log_str)
      end if
    end if
  end if
end subroutine log_error

subroutine log_warning(msg)
  character(len=*), optional :: msg
  if (log_level .ge. LOG_LEVEL_WARNING) then
    if( log_proc < 0 .or. log_proc .eq. myproc ) then
      if (present(msg)) then
        write(0,'(A)') trim(msg)
      else 
        write(0,'(A)') trim(log_str)
      end if 
    end if
  end if
end subroutine log_warning

subroutine log_info(msg)
  character(len=*), optional :: msg
  if (log_level .ge. LOG_LEVEL_INFO) then
    if( log_proc < 0 .or. log_proc .eq. myproc ) then
      if (present(msg)) then
        write(0,'(A)') trim(msg)
      else 
        write(0,'(A)') trim(log_str)
      end if
    end if
  end if
end subroutine log_info

subroutine log_debug(msg)
  character(len=*), optional :: msg
  if (log_level .ge. LOG_LEVEL_DEBUG) then
    if( log_proc < 0 .or. log_proc .eq. myproc ) then
      if (present(msg)) then
        write(0,'(A)') trim(msg)
      else 
        write(0,'(A)') trim(log_str)
      end if 
    end if
  end if
end subroutine log_debug

function str_integer(int,optional_format) result(str)
  integer, intent(in) :: int
  character(len=*), optional, intent(in) :: optional_format
  character(len=20) :: str
  if (present(optional_format)) then
    write(str,optional_format) int
  else
    write(str,*) int
  end if
end function str_integer

function str_real(re,optional_format) result(str)
  real(kind=jprb), intent(in) :: re
  character(len=*), optional :: optional_format
  character(len=20) :: str
  if (present(optional_format)) then
    write(str,optional_format) re
  else
    write(str,*) re
  end if
end function str_real

end module common_module