! (C) Copyright 2013-2015 ECMWF.

function atlas_LogChannel__ctor( cat ) result(channel)
  type(atlas_LogChannel) :: channel
  integer(c_int) :: cat
  call channel%reset_c_ptr( atlas__LogChannel_cat(cat) )
  channel%cat = cat
end function


function atlas_Logger__ctor() result(logger)
  type(atlas_Logger) :: logger
  logger%channel_error   = atlas_LogChannel(ATLAS_LOG_CATEGORY_ERROR)
  logger%channel_warning = atlas_LogChannel(ATLAS_LOG_CATEGORY_WARNING)
  logger%channel_info    = atlas_LogChannel(ATLAS_LOG_CATEGORY_INFO)
  logger%channel_debug   = atlas_LogChannel(ATLAS_LOG_CATEGORY_DEBUG)
  logger%channel_stats   = atlas_LogChannel(ATLAS_LOG_CATEGORY_STATS)
end function

subroutine LogChannel__log(this,msg,lvl,endl,flush)
  CLASS(atlas_LogChannel), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0 ; if( present(lvl)  ) opt_lvl = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_cat(this%cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_cat(this%cat,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine LogChannel__connect_stdout(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__connect_stdout(this%c_ptr())
end subroutine

subroutine LogChannel__disconnect_stdout(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__disconnect_stdout(this%c_ptr())
end subroutine

subroutine LogChannel__connect_stderr(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__connect_stderr(this%c_ptr())
end subroutine

subroutine LogChannel__disconnect_stderr(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__disconnect_stderr(this%c_ptr())
end subroutine

subroutine LogChannel__connect_fortran_unit(this,unit)
  CLASS(atlas_LogChannel) :: this
  integer, intent(in) :: unit
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__connect_fortran_unit(this%c_ptr(),unit)
end subroutine

subroutine LogChannel__disconnect_fortran_unit(this,unit)
  CLASS(atlas_LogChannel) :: this
  integer, intent(in) :: unit
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__disconnect_fortran_unit(this%c_ptr(),unit)
end subroutine

subroutine LogChannel__set_prefix(this,prefix)
  CLASS(atlas_LogChannel) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__set_prefix(this%c_ptr(),c_str_no_trim(prefix))
end subroutine

subroutine LogChannel__set_prefix_stdout(this,prefix)
  CLASS(atlas_LogChannel) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  if( .not. this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__set_prefix_stdout(this%c_ptr(),c_str_no_trim(prefix))
end subroutine

subroutine LogChannel__set_prefix_stderr(this,prefix)
  CLASS(atlas_LogChannel) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__set_prefix_stderr(this%c_ptr(),c_str_no_trim(prefix))
end subroutine

subroutine LogChannel__set_prefix_fortran_unit(this,unit,prefix)
  CLASS(atlas_LogChannel) :: this
  integer, intent(in) :: unit
  character(kind=c_char,len=*), intent(in):: prefix
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__set_prefix_fortran_unit(this%c_ptr(),unit,c_str_no_trim(prefix))
end subroutine

subroutine LogChannel__indent(this,indent)
  CLASS(atlas_LogChannel) :: this
  character(kind=c_char,len=*), intent(in), optional :: indent
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  if( present(indent) ) then
    call atlas__LogChannel__indent(this%c_ptr(),c_str_no_trim(indent))
  else
    call atlas__LogChannel__indent(this%c_ptr(),c_str_no_trim(default_indent))
  endif
end subroutine

subroutine LogChannel__indent_stdout(this,indent)
  CLASS(atlas_LogChannel) :: this
  character(kind=c_char,len=*), intent(in), optional :: indent
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  if( present(indent) ) then
    call atlas__LogChannel__indent_stdout(this%c_ptr(),c_str_no_trim(indent))
  else
    call atlas__LogChannel__indent_stdout(this%c_ptr(),c_str_no_trim(default_indent))
  endif
end subroutine

subroutine LogChannel__indent_stderr(this,indent)
  CLASS(atlas_LogChannel) :: this
  character(kind=c_char,len=*), intent(in), optional :: indent
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  if( present(indent) ) then
    call atlas__LogChannel__indent_stderr(this%c_ptr(),c_str_no_trim(indent))
  else
    call atlas__LogChannel__indent_stderr(this%c_ptr(),c_str_no_trim(default_indent))
  endif
end subroutine

subroutine LogChannel__indent_fortran_unit(this,unit,indent)
  CLASS(atlas_LogChannel) :: this
  integer, intent(in) :: unit
  character(kind=c_char,len=*), intent(in), optional :: indent
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  if( present(indent) ) then
    call atlas__LogChannel__indent_fortran_unit(this%c_ptr(),unit,c_str_no_trim(indent))
  else
    call atlas__LogChannel__indent_fortran_unit(this%c_ptr(),unit,c_str_no_trim(default_indent))
  endif
end subroutine

subroutine LogChannel__dedent(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__dedent(this%c_ptr())
end subroutine

subroutine LogChannel__dedent_stdout(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__dedent_stdout(this%c_ptr())
end subroutine

subroutine LogChannel__dedent_stderr(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__dedent_stderr(this%c_ptr())
end subroutine

subroutine LogChannel__dedent_fortran_unit(this,unit)
  CLASS(atlas_LogChannel) :: this
  integer, intent(in) :: unit
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__dedent_fortran_unit(this%c_ptr(),unit)
end subroutine

subroutine LogChannel__clear_indentation(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__clear_indentation(this%c_ptr())
end subroutine

subroutine LogChannel__clear_indentation_stdout(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__clear_indentation_stdout(this%c_ptr())
end subroutine

subroutine LogChannel__clear_indentation_stderr(this)
  CLASS(atlas_LogChannel) :: this
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__clear_indentation_stderr(this%c_ptr())
end subroutine

subroutine LogChannel__clear_indentation_fortran_unit(this,unit)
  CLASS(atlas_LogChannel) :: this
  integer, intent(in) :: unit
  if( this%is_null() ) call this%reset_c_ptr( atlas__LogChannel_cat(this%cat) )
  call atlas__LogChannel__clear_indentation_fortran_unit(this%c_ptr(),unit)
end subroutine

subroutine atlas_LogChannel__delete(this)
  class(atlas_LogChannel), intent(inout) :: this
end subroutine


subroutine atlas_LogChannel__copy(this,obj_in)
  class(atlas_LogChannel), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine


function Logger__channel(cat)
  type(atlas_LogChannel) :: Logger__channel
  integer, intent(in) :: cat
  call Logger__channel%reset_c_ptr( atlas__LogChannel_cat(cat) )
  Logger__channel%cat = cat
end function

subroutine Logger__cat(this,cat,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: cat
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 1 ; if( present(lvl)  ) opt_lvl = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_cat(cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_cat(cat,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__error(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 1 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_cat(ATLAS_LOG_CATEGORY_ERROR,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_cat(ATLAS_LOG_CATEGORY_ERROR,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__warning(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 1 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_cat(ATLAS_LOG_CATEGORY_WARNING,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_cat(ATLAS_LOG_CATEGORY_WARNING,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__info(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 1 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_cat(ATLAS_LOG_CATEGORY_INFO,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_cat(ATLAS_LOG_CATEGORY_INFO,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__set_debug(this,level)
  CLASS(atlas_Logger), intent(in) :: this
  integer , intent(in) :: level
  call atlas__log_set_debug(level)
end subroutine


subroutine Logger__debug(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 1 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_debug(opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_debug(opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__stats(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 1 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log_cat(ATLAS_LOG_CATEGORY_STATS,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log_cat(ATLAS_LOG_CATEGORY_STATS,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine


subroutine Logger__panic(this,msg)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  if( present(msg) ) then
    write(0,*) msg
  else
    write(0,*) this%msg
  endif
end subroutine

subroutine Logger__connect_fortran_unit(this,unit)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: unit
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__connect_fortran_unit(jcat,unit)
  enddo
end subroutine

subroutine Logger__disconnect_fortran_unit(this,unit)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: unit
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__disconnect_fortran_unit(jcat,unit)
  enddo
end subroutine

subroutine Logger__connect_stdout(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__connect_stdout(jcat)
  enddo
end subroutine

subroutine Logger__disconnect_stdout(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__disconnect_stdout(jcat)
  enddo
end subroutine

subroutine Logger__connect_stderr(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__connect_stderr(jcat)
  enddo
end subroutine

subroutine Logger__disconnect_stderr(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__disconnect_stderr(jcat)
  enddo
end subroutine


subroutine Logger__set_prefix_stdout(this,prefix)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  integer :: jcat
  do jcat=0,4
      call atlas__logcat__set_prefix_stdout(jcat,c_str_no_trim(prefix))
  enddo
end subroutine


subroutine Logger__set_prefix_stderr(this,prefix)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  integer :: jcat
  do jcat=0,4
      call atlas__logcat__set_prefix_stderr(jcat,c_str_no_trim(prefix))
  enddo
end subroutine

subroutine Logger__set_prefix_fortran_unit(this,unit,prefix)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: unit
  character(kind=c_char,len=*), intent(in):: prefix
  integer :: jcat
  do jcat=0,4
      call atlas__logcat__set_prefix_fortran_unit(jcat,unit,c_str_no_trim(prefix))
  enddo
end subroutine

subroutine Logger__indent(this,indent)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: indent
  integer :: jcat
  do jcat=0,4
    if( present(indent) ) then
      call atlas__logcat__indent(jcat,c_str_no_trim(indent))
    else
      call atlas__logcat__indent(jcat,c_str_no_trim(default_indent))
    endif
  enddo
end subroutine

subroutine Logger__dedent(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__dedent(jcat)
  enddo
end subroutine


subroutine Logger__clear_indentation(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__logcat__clear_indentation(jcat)
  enddo
end subroutine

subroutine atlas_Logger__delete(this)
  class(atlas_Logger), intent(inout) :: this
end subroutine

subroutine atlas_write_to_fortran_unit(unit,msg_cptr) bind(C)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int), value, intent(in) :: unit
  type(c_ptr), value, intent(in) :: msg_cptr
  character(len=:), allocatable :: msg
  msg = c_to_f_string_cptr(msg_cptr)
  write(unit,'(A)',advance='no') msg
end subroutine
