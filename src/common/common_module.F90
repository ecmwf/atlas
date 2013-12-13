module common_module

use parallel_module, only: myproc
use, intrinsic :: iso_c_binding, only: C_FLOAT, C_DOUBLE

implicit none
save

private

! Integer kind
integer, parameter, public :: jpim = selected_int_kind(9)

! Single precision
integer, parameter, public :: jprs = C_FLOAT ! selected_real_kind(4,2)
integer, parameter, public :: jprm = C_FLOAT ! selected_real_kind(6,37)

! Double precision
integer, parameter, public :: jprb = C_DOUBLE !selected_real_kind(13,300)

! Working precision = double precision
integer, parameter, public :: jprw = jprb


integer, parameter, public :: XX = 1
integer, parameter, public :: YY = 2
integer, parameter, public :: ZZ = 3

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
public :: progress_bar
public :: plot1d

integer, public :: log_level = 3
integer, public :: log_proc = -1 ! means everyone
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
  real(kind=jprw), intent(in) :: array(:)
  real(kind=jprw) :: norm
  integer :: i
  norm = 0
  do i=1,size(array)
    norm = norm + array(i)**2
  end do
  norm = sqrt(norm/real(size(array)))
end function L2norm

function Timer_elapsed(self) result(time)
  class(Timer_type), intent(inout) :: self
  real(kind=jprw) :: time
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
  real(kind=jprw), intent(in) :: re
  character(len=*), optional :: optional_format
  character(len=20) :: str
  if (present(optional_format)) then
    write(str,optional_format) re
  else
    write(str,*) re
  end if
end function str_real

subroutine progress_bar(x,xmin,xmax)
  implicit none
  real(kind=jprw), intent(in) :: x,xmin,xmax
  integer, parameter :: divisions = 51
  real(kind=jprw) :: progress_ratio
  integer, save :: prev_progress = 0
  integer :: j
    
  progress_ratio = (x-xmin)/(xmax-xmin)
  j = int(progress_ratio*divisions)
  if (j .lt. prev_progress) then
    prev_progress = 0
  end if
  if (j .gt. prev_progress) then
    if (j .eq. 1) then
      write(0,'(A)') '0%   10   20   30   40   50   60   70   80   90   100'
      write(0,'(A)') '|----|----|----|----|----|----|----|----|----|----|'
    end if
    if (j .ne. divisions) then
      write(0,'(A,$)') '*'
      flush(0)
    else 
       write(0,'(A)') '*'
    end if
    prev_progress = j
  end if
  return
end subroutine progress_bar

subroutine plot1d(width,height,val,minval,maxval)
    integer, intent(in) :: width, height
    !real(kind=jprw), intent(in)          :: val(:)
    !real(kind=jprw), intent(in), optional :: minval, maxval
    integer, intent(in)          :: val(:)
    integer, intent(in), optional :: minval, maxval
    integer :: nx, rows, cols, j, i, irow, n
    real(kind=jprw) :: vmin, vmax, scale
    rows = height
    cols = width
    nx = size(val)

    ! set or determine scaling factor for data points
    vmin= 1.0e30
    vmax=-1.0e30
    if (present(maxval)) then
      vmax = real(maxval,jprw)
    else
      do j=1,nx
        vmax = max( vmax, real(val(j),jprw) )
      end do
    end if

    if (present(minval)) then
      vmin = real(minval,jprw)
    else
      ! find absolute maximum value for scaling
      do j=1,nx
        vmin = min( vmin, real(val(j),jprw) )
      end do
    end if
    scale = real(rows-1,jprw)/(vmax-vmin+1e-6)

    ! now plot

    write(0,'(A11)',advance='no') ' '
    do j=1,cols
      write(0,'(A1)',advance='no') '-'
    end do
    write(0,'(A)') '--'

    do j=1,rows
      write(0,'(I10,A2)',advance='no') nint(vmax-(j-1)*(vmax-vmin)/real(rows-1)),'|'
      do i=1,cols
        n = 1+nint(real(i-1,jprw)*real(nx-1,jprw)/real(cols-1+1e-6,jprw))
        irow = 1+nint((vmax-real(val(n),jprw))*scale)
        if (irow.eq.j) then
          write(0,'(A1)',advance='no') '+'
        else
          write(0,'(A1)',advance='no') ' '
        end if
      end do
      write(0,'(A1)')'|'
    end do
    write(0,'(A11)',advance='no') ' '
    do i=1,cols
        write(0,'(A1)',advance='no') '-'
    end do
    write(0,'(A)') '--'

  end subroutine plot1d

end module common_module