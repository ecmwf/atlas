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
  procedure :: setup => GatherScatter__setup
  procedure :: glb_dof => GatherScatter__glb_dof
  procedure, private :: GatherScatter__gather_int32_r1_r1
  procedure, private :: GatherScatter__gather_int32_r2_r2
  procedure, private :: GatherScatter__gather_int32_r3_r3
  procedure, private :: GatherScatter__gather_real32_r1_r1
  procedure, private :: GatherScatter__gather_real32_r2_r2
  procedure, private :: GatherScatter__gather_real32_r3_r3
  procedure, private :: GatherScatter__gather_real64_r1_r1
  procedure, private :: GatherScatter__gather_real64_r2_r2
  procedure, private :: GatherScatter__gather_real64_r3_r3
  generic :: execute => &
      & GatherScatter__gather_int32_r1_r1, &
      & GatherScatter__gather_int32_r2_r2, &
      & GatherScatter__gather_int32_r3_r3, &
      & GatherScatter__gather_real32_r1_r1, &
      & GatherScatter__gather_real32_r2_r2, &
      & GatherScatter__gather_real32_r3_r3, &
      & GatherScatter__gather_real64_r1_r1, &
      & GatherScatter__gather_real64_r2_r2, &
      & GatherScatter__gather_real64_r3_r3
  generic :: gather => &
      & GatherScatter__gather_int32_r1_r1, &
      & GatherScatter__gather_int32_r2_r2, &
      & GatherScatter__gather_int32_r3_r3, &
      & GatherScatter__gather_real32_r1_r1, &
      & GatherScatter__gather_real32_r2_r2, &
      & GatherScatter__gather_real32_r3_r3, &
      & GatherScatter__gather_real64_r1_r1, &
      & GatherScatter__gather_real64_r2_r2, &
      & GatherScatter__gather_real64_r3_r3

END TYPE Gather_type
!------------------------------------------------------------------------------
