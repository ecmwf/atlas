! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: HaloExchange_type

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
      & HaloExchange__execute_real32_r1, &
      & HaloExchange__execute_real32_r2, &
      & HaloExchange__execute_real32_r3, &
      & HaloExchange__execute_real64_r1, &
      & HaloExchange__execute_real64_r2, &
      & HaloExchange__execute_real64_r3, &
      & HaloExchange__execute_real64_r4
END TYPE HaloExchange_type
!------------------------------------------------------------------------------
