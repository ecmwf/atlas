! (C) Copyright 2013-2014 ECMWF.


integer, parameter, public :: ATLAS_LOG_CATEGORY_ALL = -1
integer, parameter, public :: ATLAS_LOG_CATEGORY_ERROR = 0
integer, parameter, public :: ATLAS_LOG_CATEGORY_WARNING = 1
integer, parameter, public :: ATLAS_LOG_CATEGORY_INFO = 2
integer, parameter, public :: ATLAS_LOG_CATEGORY_DEBUG = 3
integer, parameter, public :: ATLAS_LOG_CATEGORY_STATS = 4

TYPE, extends(object_type) :: atlas_LogChannel
  character(len=1024), public :: msg
  integer, private :: cat
contains
  procedure, public :: log => LogChannel__log
  procedure, public :: connect_stdout          => LogChannel__connect_stdout
  procedure, public :: connect_stderr          => LogChannel__connect_stderr
  procedure, public :: connect_fortran_unit    => LogChannel__connect_fortran_unit

  procedure, public :: disconnect_stdout       => LogChannel__disconnect_stdout
  procedure, public :: disconnect_stderr       => LogChannel__disconnect_stderr
  procedure, public :: disconnect_fortran_unit => LogChannel__disconnect_fortran_unit

  procedure, public :: set_prefix_stdout       => LogChannel__set_prefix_stdout
  procedure, public :: set_prefix_stderr       => LogChannel__set_prefix_stderr
  procedure, public :: set_prefix_fortran_unit => LogChannel__set_prefix_fortran_unit

ENDTYPE

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_Logger

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
!   cat     : Log to channel by category index

! Author :
! ------
!   Sept-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------

  character(len=1024), public :: msg

  type(ATLAS_LogChannel) :: channel_error
  type(ATLAS_LogChannel) :: channel_warning
  type(ATLAS_LogChannel) :: channel_info
  type(ATLAS_LogChannel) :: channel_debug
  type(ATLAS_LogChannel) :: channel_stats

contains

  procedure, public :: init => Logger__init

  procedure, public :: channel => Logger__channel
  procedure, public :: set_debug => Logger__set_debug
  procedure, public :: debug => Logger__debug
  procedure, public :: info => Logger__info
  procedure, public :: warning => Logger__warning
  procedure, public :: error => Logger__error
  procedure, public :: panic => Logger__panic
  procedure, public :: stats => Logger__stats
  procedure, public :: cat => Logger__cat

  procedure, public :: connect_stdout          => Logger__connect_stdout
  procedure, public :: connect_stderr          => Logger__connect_stderr
  procedure, public :: connect_fortran_unit    => Logger__connect_fortran_unit
  procedure, public :: disconnect_fortran_unit => Logger__disconnect_fortran_unit
  procedure, public :: disconnect_stdout       => Logger__disconnect_stdout
  procedure, public :: disconnect_stderr       => Logger__disconnect_stderr
  procedure, public :: set_prefix_stdout       => Logger__set_prefix_stdout
  procedure, public :: set_prefix_stderr       => Logger__set_prefix_stderr
  procedure, public :: set_prefix_fortran_unit => Logger__set_prefix_fortran_unit

END TYPE

!------------------------------------------------------------------------------


