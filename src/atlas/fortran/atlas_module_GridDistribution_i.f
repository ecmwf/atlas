! (C) Copyright 2013-2014 ECMWF.



!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_GridDistribution

! Purpose :
! -------
!   *GridDistribution* : Object passed from atlas to inspect grid distribution

! Methods :
! -------

! Author :
! ------
!   12-Mar-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_GridDistribution

!------------------------------------------------------------------------------

interface new_atlas_GridDistribution
  module procedure new_atlas_GridDistribution
end interface
