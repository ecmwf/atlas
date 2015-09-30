! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_Checksum

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
  procedure, private :: Checksum__setup32
  procedure, private :: Checksum__setup64
  procedure, private :: Checksum__execute_int32_r1
  procedure, private :: Checksum__execute_int32_r2
  procedure, private :: Checksum__execute_int32_r3
  procedure, private :: Checksum__execute_real32_r1
  procedure, private :: Checksum__execute_real32_r2
  procedure, private :: Checksum__execute_real32_r3
  procedure, private :: Checksum__execute_real64_r1
  procedure, private :: Checksum__execute_real64_r3
  procedure, private :: Checksum__execute_real64_r2
  generic :: setup => &
      & Checksum__setup32, &
      & Checksum__setup64
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

END TYPE atlas_Checksum

!------------------------------------------------------------------------------

interface atlas_Checksum
  module procedure atlas_Checksum__ctor
end interface

!------------------------------------------------------------------------------
