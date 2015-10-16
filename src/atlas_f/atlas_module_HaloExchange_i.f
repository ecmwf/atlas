! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_HaloExchange

! Purpose :
! -------
!   *HaloExchange* :

! Methods :
! -------
!   setup : Setup using arrays detailing proc and glb_idx, bounds and parbound
!   execute : Do the halo exchange

! Author :
! ------
!   17-Dec-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: setup => HaloExchange__setup
  procedure, private :: HaloExchange__execute_int32_r1
  procedure, private :: HaloExchange__execute_int32_r2
  procedure, private :: HaloExchange__execute_int32_r3
  procedure, private :: HaloExchange__execute_int64_r1
  procedure, private :: HaloExchange__execute_int64_r2
  procedure, private :: HaloExchange__execute_int64_r3
  procedure, private :: HaloExchange__execute_real32_r1
  procedure, private :: HaloExchange__execute_real32_r2
  procedure, private :: HaloExchange__execute_real32_r3
  procedure, private :: HaloExchange__execute_real64_r1
  procedure, private :: HaloExchange__execute_real64_r2
  procedure, private :: HaloExchange__execute_real64_r3
  procedure, private :: HaloExchange__execute_real64_r4
  generic :: execute => &
      & HaloExchange__execute_int32_r1, &
      & HaloExchange__execute_int32_r2, &
      & HaloExchange__execute_int32_r3, &
      & HaloExchange__execute_int64_r1, &
      & HaloExchange__execute_int64_r2, &
      & HaloExchange__execute_int64_r3, &
      & HaloExchange__execute_real32_r1, &
      & HaloExchange__execute_real32_r2, &
      & HaloExchange__execute_real32_r3, &
      & HaloExchange__execute_real64_r1, &
      & HaloExchange__execute_real64_r2, &
      & HaloExchange__execute_real64_r3, &
      & HaloExchange__execute_real64_r4

  procedure, public :: delete => atlas_HaloExchange__delete

END TYPE atlas_HaloExchange
!------------------------------------------------------------------------------

interface atlas_HaloExchange
  module procedure atlas_HaloExchange__ctor
end interface

!------------------------------------------------------------------------------

