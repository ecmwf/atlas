
module atlas_Error_module

implicit none
private

public :: atlas_CodeLocation
public :: atlas_code_location_str
public :: atlas_code_location
public :: atlas_abort
public :: atlas_throw_exception
public :: atlas_throw_notimplemented
public :: atlas_throw_outofrange
public :: atlas_throw_seriousbug
public :: atlas_throw_usererror
public :: atlas_throw_assertionfailed
public :: atlas_err_code
public :: atlas_err_msg
public :: atlas_err_set_aborts
public :: atlas_err_set_throws
public :: atlas_err_set_backtrace
public :: atlas_err
public :: atlas_noerr
public :: atlas_err_clear
public :: atlas_err_success

! Error codes
integer, parameter, public ::      &
  atlas_err_cleared         =  1,  &
  atlas_err_noerr           =  0,   &
  atlas_err_exception       = -1,   &
  atlas_err_usererror       = -2,   &
  atlas_err_seriousbug      = -3,   &
  atlas_err_notimplemented  = -4,   &
  atlas_err_assertionfailed = -5,   &
  atlas_err_badparameter    = -6,   &
  atlas_err_outofrange      = -7,   &
  atlas_err_stop            = -100, &
  atlas_err_abort           = -101, &
  atlas_err_cancel          = -102, &
  atlas_err_readerror       = -200, &
  atlas_err_writeerror      = -201, &
  atlas_err_unknown         = -999

integer, private, parameter :: ATLAS_CODELOCATION_FILE_STRLEN     = 1024
integer, private, parameter :: ATLAS_CODELOCATION_FUNCTION_STRLEN = 1024

TYPE :: atlas_CodeLocation
  integer :: line
  character(len=ATLAS_CODELOCATION_FILE_STRLEN) :: file
  character(len=ATLAS_CODELOCATION_FILE_STRLEN) :: function
contains
  procedure :: str => CodeLocation__str
ENDTYPE

interface atlas_code_location_str
  module procedure code_location_str_FILE_LINE
end interface

interface atlas_abort
  module procedure atlas_abort_null
  module procedure atlas_abort_msg
  module procedure atlas_abort_msg_loc
end interface atlas_abort

interface atlas_throw_exception
  module procedure atlas_throw_exception_msg
  module procedure atlas_throw_exception_msg_loc
  module procedure atlas_throw_exception_loc
end interface atlas_throw_exception

interface atlas_throw_notimplemented
  module procedure atlas_throw_notimplemented_msg
  module procedure atlas_throw_notimplemented_loc
  module procedure atlas_throw_notimplemented_msg_loc
end interface atlas_throw_notimplemented

interface atlas_throw_outofrange
  module procedure atlas_throw_outofrange_msg
  module procedure atlas_throw_outofrange_loc
  module procedure atlas_throw_outofrange_msg_loc
  module procedure atlas_throw_outofrange_range
  module procedure atlas_throw_outofrange_range_loc
end interface atlas_throw_outofrange

interface atlas_throw_seriousbug
  module procedure atlas_throw_seriousbug_msg
  module procedure atlas_throw_seriousbug_loc
  module procedure atlas_throw_seriousbug_msg_loc
end interface atlas_throw_seriousbug

interface atlas_throw_usererror
  module procedure atlas_throw_usererror_msg
  module procedure atlas_throw_usererror_loc
  module procedure atlas_throw_usererror_msg_loc
end interface atlas_throw_usererror

interface atlas_throw_assertionfailed
  module procedure atlas_throw_assertionfailed_msg
  module procedure atlas_throw_assertionfailed_loc
  module procedure atlas_throw_assertionfailed_msg_loc
end interface atlas_throw_assertionfailed


interface atlas_code_location
  module procedure code_location_null
  module procedure code_location_file_line
  module procedure code_location_file_line_func
end interface atlas_code_location

!------------------------------------------------------------------------------
!========================================================
contains
!========================================================



function atlas_err_code()
  use atlas_errorhandling_c_binding
  use, intrinsic :: iso_c_binding, only: c_ptr
  integer :: atlas_err_code
  atlas_err_code = atlas__Error_code()
end function

function atlas_err_msg()
  use, intrinsic :: iso_c_binding, only: c_ptr
  use atlas_errorhandling_c_binding
  use fckit_c_interop_module , only:  c_ptr_to_string
  type(c_ptr) :: msg_cptr
  character(len=:), allocatable :: atlas_err_msg
  msg_cptr = atlas__Error_msg()
  atlas_err_msg = c_ptr_to_string(msg_cptr)
end function

subroutine atlas_err_set_aborts( aborts )
  use atlas_errorhandling_c_binding
  logical, intent(in) :: aborts
  if( aborts ) then
    call atlas__Error_set_aborts(1)
  else
    call atlas__Error_set_aborts(0)
  endif
end subroutine

subroutine atlas_err_set_throws( throws )
  use atlas_errorhandling_c_binding
  logical, intent(in) :: throws
  if( throws ) then
    call atlas__Error_set_throws(1)
  else
    call atlas__Error_set_throws(0)
  endif
end subroutine

subroutine atlas_err_set_backtrace( backtrace )
  use atlas_errorhandling_c_binding
  logical, intent(in) :: backtrace
  if( backtrace ) then
    call atlas__Error_set_backtrace(1)
  else
    call atlas__Error_set_backtrace(0)
  endif
end subroutine


function CodeLocation__str(self) result( str )
  class(atlas_CodeLocation) :: self
  character(len(self%file)+5) :: str
  write(str,'(A,A1,I4)') self%file,":",self%line
end function

function code_location_str_FILE_LINE(file,line) result( str )
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  character(len(file)+5) :: str
  type(atlas_CodeLocation) :: code_location
  code_location = code_location_file_line(file,line)
  str = code_location%str()
end function


function code_location_null() result( code_location )
  type(atlas_CodeLocation) :: code_location
  code_location%file = ""
  code_location%line = 0
  code_location%function = ""
end function

function code_location_file_line(file,line) result( code_location )
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  type(atlas_CodeLocation) :: code_location
  code_location%file = file
  code_location%line = line
  code_location%function = ""
end function

function code_location_file_line_func(file,line,func) result( code_location )
  character(len=*), intent(in) :: file, func
  integer, intent(in) :: line
  type(atlas_CodeLocation) :: code_location
  code_location%file = file
  code_location%line = line
  code_location%function = func
end function

subroutine atlas_abort_null()
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  call atlas__abort(c_str(""),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_abort_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__abort(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_abort_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__abort(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_exception_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_exception(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_exception_loc(code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_exception(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_exception_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_exception(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_notimplemented_loc(code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_notimplemented(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_notimplemented_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_notimplemented(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_notimplemented_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_notimplemented(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_outofrange_loc(code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_outofrange(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_outofrange_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_outofrange(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_outofrange_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_outofrange(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_outofrange_range_loc(arrayname,idx,max,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: arrayname
  integer, intent(in) :: idx, max
  type(atlas_CodeLocation), intent(in) :: code_loc
  character(len=80) :: msg
  write(msg,'(A,I0,A,I0,A,A)') "Index ",idx," is greater than maximum ",max," in array ",arrayname
  call atlas__throw_outofrange(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_outofrange_range(arrayname,idx,max)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: arrayname
  integer, intent(in) :: idx, max
  character(len=80) :: msg
  write(msg,'(A,I0,A,I0,A,A)') "Index ",idx," is greater than maximum ",max," in array ",arrayname
  call atlas__throw_outofrange(c_str(msg),c_str(""),0,c_str(""))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_usererror_loc(code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_usererror(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_usererror_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_usererror(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_usererror_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_usererror(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_assertionfailed_loc(code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_assertionfailed(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_assertionfailed_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_assertionfailed(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_assertionfailed_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_assertionfailed(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_seriousbug_loc(code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_seriousbug(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_seriousbug_msg(msg)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_seriousbug(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_seriousbug_msg_loc(msg,code_loc)
  use fckit_c_interop_module, only:  c_str
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(atlas_CodeLocation), intent(in) :: code_loc
  call atlas__throw_seriousbug(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------
subroutine atlas_err_clear()
  use atlas_errorhandling_c_binding
  call atlas__Error_clear()
end subroutine

subroutine atlas_err_success()
  use atlas_errorhandling_c_binding
  call atlas__Error_success()
end subroutine

function atlas_noerr()
  logical :: atlas_noerr
  if( atlas_err_code() == atlas_err_noerr ) then
    atlas_noerr = .True.
  else
    atlas_noerr = .False.
  endif
end function

function atlas_err()
  logical :: atlas_err
  if( atlas_err_code() /= atlas_err_noerr ) then
    atlas_err = .True.
  else
    atlas_err = .False.
  endif
end function

end module atlas_Error_module

