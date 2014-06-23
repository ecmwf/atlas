! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: Gather_type

! Purpose :
! -------
!   *Gather* :

! Methods :
! -------
!   setup : Setup using arrays detailing proc, glb_idx, remote_idx, max_glb_idx
!   execute : Do the gather

! Author :
! ------
!   17-Dec-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: setup => Gather__setup
  procedure :: glb_dof => Gather__glb_dof
  procedure, private :: Gather__execute_int32_r1_r1
  procedure, private :: Gather__execute_int32_r2_r2
  procedure, private :: Gather__execute_int32_r3_r3
  procedure, private :: Gather__execute_real32_r1_r1
  procedure, private :: Gather__execute_real32_r2_r2
  procedure, private :: Gather__execute_real32_r3_r3
  procedure, private :: Gather__execute_real64_r1_r1
  procedure, private :: Gather__execute_real64_r2_r2
  procedure, private :: Gather__execute_real64_r3_r3
  generic :: execute => &
      & Gather__execute_int32_r1_r1, &
      & Gather__execute_int32_r2_r2, &
      & Gather__execute_int32_r3_r3, &
      & Gather__execute_real32_r1_r1, &
      & Gather__execute_real32_r2_r2, &
      & Gather__execute_real32_r3_r3, &
      & Gather__execute_real64_r1_r1, &
      & Gather__execute_real64_r2_r2, &
      & Gather__execute_real64_r3_r3
END TYPE Gather_type
!------------------------------------------------------------------------------
