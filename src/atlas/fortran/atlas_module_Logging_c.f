! (C) Copyright 2013-2014 ECMWF.

subroutine Logger__cat(this,cat,msg,lvl,endl,flush)
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
    call atlas__log(cat,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__error(this,msg,lvl,endl,flush)
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0 ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(this%cat_error,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(this%cat_error,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__warning(this,msg,lvl,endl,flush)
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(this%cat_warning,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(this%cat_warning,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__info(this,msg,lvl,endl,flush)
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(this%cat_info,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(this%cat_info,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
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
    call atlas__log_debug(opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine

subroutine Logger__stats(this,msg,lvl,endl,flush)
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  integer, intent(in), optional :: lvl
  logical, intent(in), optional :: endl, flush
  integer :: opt_lvl, opt_endl, opt_flush
  opt_lvl   = 0      ; if( present(lvl)  ) opt_lvl   = lvl
  opt_endl  = 1 ; if( present(endl) ) then; if( .not. endl )  opt_endl  = 0; endif
  opt_flush = 1 ; if( present(flush)) then; if( .not. flush ) opt_flush = 0; endif
  if (present(msg)) then
    call atlas__log(this%cat_stats,opt_lvl,c_str(trim(msg)),opt_endl,opt_flush)
  else
    call atlas__log(this%cat_stats,opt_lvl,c_str(trim(this%msg)),opt_endl,opt_flush)
  end if
end subroutine


subroutine Logger__panic(this,msg)
  CLASS(Logger_t), intent(in) :: this
  character(kind=c_char,len=*), intent(in), optional :: msg
  if( present(msg) ) then
    write(0,*) msg
  else
    write(0,*) this%msg
  endif
end subroutine

subroutine Logger__connect_fortran_unit(this,cat,unit)
  CLASS(Logger_t), intent(in) :: this
  integer, intent(in) :: cat
  integer, intent(in) :: unit
  integer :: jcat
  if( cat == ATLAS_LOG_CAT_ALL ) then
    do jcat=0,4
      call atlas__cat__connect_fortran_unit(jcat,unit)
    enddo
  else
    call atlas__cat__connect_fortran_unit(cat,unit)
  endif
end subroutine


subroutine atlas_write_to_fortran_unit(unit,msg_cptr) bind(C)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int), value, intent(in) :: unit
  type(c_ptr), value, intent(in) :: msg_cptr
  character(len=:), allocatable :: msg
  msg = c_to_f_string_cptr(msg_cptr)
  write(unit,'(A)',advance='no') msg
end subroutine
