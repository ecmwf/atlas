! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: Logger_t

! Purpose :
! -------
!   *Logger* :

! Methods :
! -------
!   debug   : Log to debug channel
!   info    : Log to info channel
!   warning : Log to warning channel
!   error   : Log to error channel
!   stats   : Log to stats channel
!   panic   : Log to panic channel

! Author :
! ------
!   Sept-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------

  character(len=1024), public :: message

contains

  procedure, public :: set_debug => Logger__set_debug
  procedure, public :: debug => Logger__debug
  procedure, public :: info => Logger__info
  procedure, public :: warning => Logger__warning
  procedure, public :: error => Logger__error
  procedure, public :: panic => Logger__panic
  procedure, public :: stats => Logger__stats
  procedure, public :: log => Logger__log
END TYPE
!------------------------------------------------------------------------------
