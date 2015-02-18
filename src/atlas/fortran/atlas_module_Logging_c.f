! (C) Copyright 2013-2014 ECMWF.

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
    call atlas__log(this%cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(this%cat,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine LogChannel__connect_stdout(this)
  CLASS(atlas_LogChannel), intent(in) :: this
  call atlas__Channel__connect_stdout(this%private%object)
end subroutine

subroutine LogChannel__disconnect_stdout(this)
  CLASS(atlas_LogChannel), intent(in) :: this
  call atlas__Channel__disconnect_stdout(this%private%object)
end subroutine

subroutine LogChannel__connect_stderr(this)
  CLASS(atlas_LogChannel), intent(in) :: this
  call atlas__Channel__connect_stderr(this%private%object)
end subroutine

subroutine LogChannel__disconnect_stderr(this)
  CLASS(atlas_LogChannel), intent(in) :: this
  call atlas__Channel__disconnect_stderr(this%private%object)
end subroutine

subroutine LogChannel__connect_fortran_unit(this,unit)
  CLASS(atlas_LogChannel), intent(in) :: this
  integer, intent(in) :: unit
  call atlas__Channel__connect_fortran_unit(this%private%object,unit)
end subroutine

subroutine LogChannel__disconnect_fortran_unit(this,unit)
  CLASS(atlas_LogChannel), intent(in) :: this
  integer, intent(in) :: unit
  call atlas__Channel__disconnect_fortran_unit(this%private%object,unit)
end subroutine


subroutine LogChannel__set_prefix_stdout(this,prefix)
  CLASS(atlas_LogChannel), intent(in) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  call atlas__Channel__set_prefix_stdout(this%private%object,c_str_no_trim(prefix))
end subroutine

subroutine LogChannel__set_prefix_stderr(this,prefix)
  CLASS(atlas_LogChannel), intent(in) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  call atlas__Channel__set_prefix_stderr(this%private%object,c_str_no_trim(prefix))
end subroutine

subroutine LogChannel__set_prefix_fortran_unit(this,unit,prefix)
  CLASS(atlas_LogChannel), intent(in) :: this
  integer, intent(in) :: unit
  character(kind=c_char,len=*), intent(in):: prefix
  call atlas__Channel__set_prefix_fortran_unit(this%private%object,unit,c_str_no_trim(prefix))
end subroutine

function Logger__channel(this,cat)
  type(atlas_LogChannel) :: Logger__channel
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: cat
  Logger__channel%private%object = atlas__log_channel(cat)
  Logger__channel%cat = cat
end function

subroutine Logger__cat(this,cat,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: cat
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0 ; if( present(lvl)  ) opt_lvl = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(cat,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__error(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(ATLAS_LOG_CATEGORY_ERROR,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(ATLAS_LOG_CATEGORY_ERROR,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__warning(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(ATLAS_LOG_CATEGORY_WARNING,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(ATLAS_LOG_CATEGORY_WARNING,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__info(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(ATLAS_LOG_CATEGORY_INFO,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(ATLAS_LOG_CATEGORY_INFO,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__set_debug(this,level)
  CLASS(atlas_Logger), intent(in) :: this
  integer , intent(in) :: level
  call atlas__log_debug_set_level(level)
end subroutine


subroutine Logger__debug(this,msg,lvl,endl,flush)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
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
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(ATLAS_LOG_CATEGORY_STATS,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(ATLAS_LOG_CATEGORY_STATS,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
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
    call atlas__cat__connect_fortran_unit(jcat,unit)
  enddo
end subroutine

subroutine Logger__disconnect_fortran_unit(this,unit)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: unit
  integer :: jcat
  do jcat=0,4
    call atlas__cat__disconnect_fortran_unit(jcat,unit)
  enddo
end subroutine

subroutine Logger__connect_stdout(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__cat__connect_stdout(jcat)
  enddo
end subroutine

subroutine Logger__disconnect_stdout(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__cat__disconnect_stdout(jcat)
  enddo
end subroutine

subroutine Logger__connect_stderr(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__cat__connect_stderr(jcat)
  enddo
end subroutine

subroutine Logger__disconnect_stderr(this)
  CLASS(atlas_Logger), intent(in) :: this
  integer :: jcat
  do jcat=0,4
    call atlas__cat__disconnect_stderr(jcat)
  enddo
end subroutine


subroutine Logger__set_prefix_stdout(this,prefix)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  integer :: jcat
  do jcat=0,4
      call atlas__cat__set_prefix_stdout(jcat,c_str_no_trim(prefix))
  enddo
end subroutine


subroutine Logger__set_prefix_stderr(this,prefix)
  CLASS(atlas_Logger), intent(in) :: this
  character(kind=c_char,len=*), intent(in):: prefix
  integer :: jcat
  do jcat=0,4
      call atlas__cat__set_prefix_stderr(jcat,c_str_no_trim(prefix))
  enddo
end subroutine

subroutine Logger__set_prefix_fortran_unit(this,unit,prefix)
  CLASS(atlas_Logger), intent(in) :: this
  integer, intent(in) :: unit
  character(kind=c_char,len=*), intent(in):: prefix
  integer :: jcat
  do jcat=0,4
      call atlas__cat__set_prefix_fortran_unit(jcat,unit,c_str_no_trim(prefix))
  enddo
end subroutine

subroutine Logger__init(this)
  CLASS(atlas_Logger), intent(inout) :: this
  this%channel_error   = this%channel(ATLAS_LOG_CATEGORY_ERROR)
  this%channel_warning = this%channel(ATLAS_LOG_CATEGORY_WARNING)
  this%channel_info    = this%channel(ATLAS_LOG_CATEGORY_INFO)
  this%channel_debug   = this%channel(ATLAS_LOG_CATEGORY_DEBUG)
  this%channel_stats   = this%channel(ATLAS_LOG_CATEGORY_STATS)
end subroutine


subroutine atlas_write_to_fortran_unit(unit,msg_cptr) bind(C)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int), value, intent(in) :: unit
  type(c_ptr), value, intent(in) :: msg_cptr
  character(len=:), allocatable :: msg
  msg = c_to_f_string_cptr(msg_cptr)
  write(unit,'(A)',advance='no') msg
end subroutine
