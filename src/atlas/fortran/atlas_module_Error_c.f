! (C) Copyright 2013-2014 ECMWF.


function atlas_err_code()
  use atlas_errorhandling_c_binding
  use iso_c_binding, only: c_ptr
  integer :: atlas_err_code
  atlas_err_code = atlas__Error_code()
end function

function atlas_err_msg()
  use atlas_errorhandling_c_binding
  use iso_c_binding, only: c_ptr
  type(c_ptr) :: msg_cptr
  character(len=:), allocatable :: atlas_err_msg
  msg_cptr = atlas__Error_msg()
  atlas_err_msg = c_to_f_string_cptr(msg_cptr)
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

function code_location_null() result( code_location )
  type(CodeLocation_type) :: code_location
  code_location%file = ""
  code_location%line = 0
  code_location%function = ""
end function

function code_location_file_line(file,line) result( code_location )
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  type(CodeLocation_type) :: code_location
  code_location%file = file
  code_location%line = line
  code_location%function = ""
end function

function code_location_file_line_func(file,line,func) result( code_location )
  character(len=*), intent(in) :: file, func
  integer, intent(in) :: line
  type(CodeLocation_type) :: code_location
  code_location%file = file
  code_location%line = line
  code_location%function = func
end function

subroutine atlas_abort_null()
  use atlas_errorhandling_c_binding
  call atlas__abort(c_str(""),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_abort_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__abort(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_abort_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__abort(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_exception_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_exception(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_exception_loc(code_loc)
  use atlas_errorhandling_c_binding
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_exception(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_exception_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_exception(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_notimplemented_loc(code_loc)
  use atlas_errorhandling_c_binding
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_notimplemented(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_notimplemented_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_notimplemented(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_notimplemented_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_notimplemented(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_outofrange_loc(code_loc)
  use atlas_errorhandling_c_binding
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_outofrange(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_outofrange_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_outofrange(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_outofrange_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_outofrange(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_outofrange_range_loc(arrayname,idx,max,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: arrayname
  integer, intent(in) :: idx, max
  type(CodeLocation_type), intent(in) :: code_loc
  character(len=80) :: msg
  write(msg,'(A,I0,A,I0,A,A)') "Index ",idx," is greater than maximum ",max," in array ",arrayname
  call atlas__throw_outofrange(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_outofrange_range(arrayname,idx,max)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: arrayname
  integer, intent(in) :: idx, max
  character(len=80) :: msg
  write(msg,'(A,I0,A,I0,A,A)') "Index ",idx," is greater than maximum ",max," in array ",arrayname
  call atlas__throw_outofrange(c_str(msg),c_str(""),0,c_str(""))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_usererror_loc(code_loc)
  use atlas_errorhandling_c_binding
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_usererror(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_usererror_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_usererror(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_usererror_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_usererror(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_assertionfailed_loc(code_loc)
  use atlas_errorhandling_c_binding
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_assertionfailed(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_assertionfailed_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_assertionfailed(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_assertionfailed_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_assertionfailed(c_str(msg),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

! -----------------------------------------------------------------------------

subroutine atlas_throw_seriousbug_loc(code_loc)
  use atlas_errorhandling_c_binding
  type(CodeLocation_type), intent(in) :: code_loc
  call atlas__throw_seriousbug(c_str(""),c_str(code_loc%file),code_loc%line,c_str(code_loc%function))
end subroutine

subroutine atlas_throw_seriousbug_msg(msg)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  call atlas__throw_seriousbug(c_str(msg),c_str(""),0,c_str(""))
end subroutine

subroutine atlas_throw_seriousbug_msg_loc(msg,code_loc)
  use atlas_errorhandling_c_binding
  character(len=*), intent(in) :: msg
  type(CodeLocation_type), intent(in) :: code_loc
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


subroutine atlas_error_example()
  use atlas_errorhandling_c_binding
  call atlas__error_example()
end subroutine
