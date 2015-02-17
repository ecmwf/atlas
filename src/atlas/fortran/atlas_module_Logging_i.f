! (C) Copyright 2013-2014 ECMWF.


integer, parameter, private :: ATLAS_LOG_CAT_ALL = -1
integer, parameter, private :: ATLAS_LOG_CAT_ERROR = 0
integer, parameter, private :: ATLAS_LOG_CAT_WARNING = 1
integer, parameter, private :: ATLAS_LOG_CAT_INFO = 2
integer, parameter, private :: ATLAS_LOG_CAT_DEBUG = 3
integer, parameter, private :: ATLAS_LOG_CAT_STATS = 4

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
!   cat     : Log to channel by category index

! Author :
! ------
!   Sept-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------

  character(len=1024), public :: msg

  integer :: cat_all      = -1
  integer :: cat_error    =  0
  integer :: cat_warning  =  1
  integer :: cat_info     =  2
  integer :: cat_debug    =  3
  integer :: cat_stats    =  4

contains

  procedure, public :: set_debug => Logger__set_debug
  procedure, public :: debug => Logger__debug
  procedure, public :: info => Logger__info
  procedure, public :: warning => Logger__warning
  procedure, public :: error => Logger__error
  procedure, public :: panic => Logger__panic
  procedure, public :: stats => Logger__stats
  procedure, public :: cat => Logger__cat
  procedure, public :: connect_fortran_unit => Logger__connect_fortran_unit
END TYPE

!------------------------------------------------------------------------------
