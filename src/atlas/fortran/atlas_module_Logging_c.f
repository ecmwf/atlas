! (C) Copyright 2013-2014 ECMWF.

subroutine Logger__log(this,cat,msg,lvl,endl,flush)
  CLASS(Logger_t), intent(in) :: this
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
    call atlas__log(cat,opt_lvl,c_str(trim(this%message)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__error(this,msg,lvl,endl,flush)
  integer, parameter :: cat = 0
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(cat,opt_lvl,c_str(trim(this%message)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__warning(this,msg,lvl,endl,flush)
  integer, parameter :: cat = 1
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(cat,opt_lvl,c_str(trim(this%message)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__info(this,msg,lvl,endl,flush)
  integer, parameter :: cat = 2
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(cat,opt_lvl,c_str(trim(this%message)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__set_debug(this,level)
  CLASS(Logger_t), intent(in) :: this
  integer , intent(in) :: level
  call atlas__log_debug_set_level(level)
end subroutine


subroutine Logger__debug(this,msg,lvl,endl,flush)
  CLASS(Logger_t), intent(in) :: this
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
    call atlas__log_debug(opt_lvl,c_str(trim(this%message)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__stats(this,msg,lvl,endl,flush)
  integer, parameter :: cat = 4
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(cat,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(cat,opt_lvl,c_str(trim(this%message)),opt_endl,opt_flush)
  end if
end subroutine


subroutine Logger__panic(this,msg)
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  if( present(msg) ) then
    write(0,*) msg
  else
    write(0,*) this%message
  endif
end subroutine
