! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: Checksum_type

! Purpose :
! -------
!   *Checksum* :

! Methods :
! -------
!   setup : Setup using arrays detailing proc, glb_idx, remote_idx, max_glb_idx
!   execute : Do the Checksum

! Author :
! ------
!   27-Jun-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: setup => Checksum__setup
  procedure, private :: Checksum__execute_int32_r1
  procedure, private :: Checksum__execute_int32_r2
  procedure, private :: Checksum__execute_int32_r3
  procedure, private :: Checksum__execute_real32_r1
  procedure, private :: Checksum__execute_real32_r2
  procedure, private :: Checksum__execute_real32_r3
  procedure, private :: Checksum__execute_real64_r1
  procedure, private :: Checksum__execute_real64_r3
  procedure, private :: Checksum__execute_real64_r2
  generic :: execute => &
      & Checksum__execute_int32_r1, &
      & Checksum__execute_int32_r2, &
      & Checksum__execute_int32_r3, &
      & Checksum__execute_real32_r1, &
      & Checksum__execute_real32_r2, &
      & Checksum__execute_real32_r3, &
      & Checksum__execute_real64_r1, &
      & Checksum__execute_real64_r2, &
      & Checksum__execute_real64_r3

END TYPE Checksum_type
!------------------------------------------------------------------------------
