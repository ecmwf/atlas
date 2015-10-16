! (C) Copyright 2013-2015 ECMWF.


integer, parameter, public :: ATLAS_LOG_CATEGORY_ALL     = -1
integer, parameter, public :: ATLAS_LOG_CATEGORY_ERROR   = 0
integer, parameter, public :: ATLAS_LOG_CATEGORY_WARNING = 1
integer, parameter, public :: ATLAS_LOG_CATEGORY_INFO    = 2
integer, parameter, public :: ATLAS_LOG_CATEGORY_DEBUG   = 3
integer, parameter, public :: ATLAS_LOG_CATEGORY_STATS   = 4

character(len=4), parameter, private :: default_indent = "    "


TYPE, extends(atlas_object) :: atlas_LogChannel
  character(len=1024), public :: msg = ""
  integer, private :: cat
contains
  procedure, public :: log => LogChannel__log
  procedure, public :: connect_stdout          => LogChannel__connect_stdout
  procedure, public :: connect_stderr          => LogChannel__connect_stderr
  procedure, public :: connect_fortran_unit    => LogChannel__connect_fortran_unit

  procedure, public :: disconnect_stdout       => LogChannel__disconnect_stdout
  procedure, public :: disconnect_stderr       => LogChannel__disconnect_stderr
  procedure, public :: disconnect_fortran_unit => LogChannel__disconnect_fortran_unit

  procedure, public :: set_prefix              => LogChannel__set_prefix
  procedure, public :: set_prefix_stdout       => LogChannel__set_prefix_stdout
  procedure, public :: set_prefix_stderr       => LogChannel__set_prefix_stderr
  procedure, public :: set_prefix_fortran_unit => LogChannel__set_prefix_fortran_unit

  procedure, public :: indent                  => LogChannel__indent
  procedure, public :: indent_stdout           => LogChannel__indent_stdout
  procedure, public :: indent_stderr           => LogChannel__indent_stderr
  procedure, public :: indent_fortran_unit     => LogChannel__indent_fortran_unit

  procedure, public :: dedent                  => LogChannel__dedent
  procedure, public :: dedent_stdout           => LogChannel__dedent_stdout
  procedure, public :: dedent_stderr           => LogChannel__dedent_stderr
  procedure, public :: dedent_fortran_unit     => LogChannel__dedent_fortran_unit

  procedure, public :: clear_indentation              => LogChannel__clear_indentation
  procedure, public :: clear_indentation_stdout       => LogChannel__clear_indentation_stdout
  procedure, public :: clear_indentation_stderr       => LogChannel__clear_indentation_stderr
  procedure, public :: clear_indentation_fortran_unit => LogChannel__clear_indentation_fortran_unit

  procedure, public :: delete => atlas_LogChannel__delete


ENDTYPE

interface atlas_LogChannel
  module procedure atlas_LogChannel__ctor
end interface


!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_Logger

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

  character(len=1024), public :: msg=""

  type(ATLAS_LogChannel) :: channel_error
  type(ATLAS_LogChannel) :: channel_warning
  type(ATLAS_LogChannel) :: channel_info
  type(ATLAS_LogChannel) :: channel_debug
  type(ATLAS_LogChannel) :: channel_stats

contains

  procedure, public, nopass :: channel => Logger__channel
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

  procedure, public :: indent            => Logger__indent
  procedure, public :: dedent            => Logger__dedent
  procedure, public :: clear_indentation => Logger__clear_indentation

  procedure, public :: delete => atlas_Logger__delete

END TYPE

interface atlas_Logger
  module procedure atlas_Logger__ctor
end interface


!------------------------------------------------------------------------------


