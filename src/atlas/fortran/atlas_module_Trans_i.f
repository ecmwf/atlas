! (C) Copyright 2013-2015 ECMWF.



!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_Trans

! Purpose :
! -------
!   *Trans* : To do spectral transforms

! Methods :
! -------

! Author :
! ------
!   12-Mar-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: handle => atlas_Trans__handle
END TYPE atlas_Trans

!------------------------------------------------------------------------------

interface new_atlas_Trans
  module procedure new_atlas_Trans
end interface
