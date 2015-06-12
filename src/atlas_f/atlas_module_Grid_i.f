! (C) Copyright 2013-2014 ECMWF.



!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_ReducedGrid

! Purpose :
! -------
!   *ReducedGrid* : Object Grid specifications for Reduced Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: npts => ReducedGrid__npts
  procedure :: nlat => ReducedGrid__nlat
  procedure :: nlon => ReducedGrid__nlon
  procedure :: lat => ReducedGrid__latitudes
END TYPE atlas_ReducedGrid

interface atlas_ReducedGrid
  module procedure atlas_ReducedGrid__ctor_id
end interface

interface atlas_GaussianGrid
  module procedure atlas_GaussianGrid__ctor
end interface

interface atlas_ReducedGaussianGrid
  module procedure atlas_ReducedGaussianGrid__ctor
end interface

interface atlas_LonLatGrid
  module procedure atlas_LonLatGrid__ctor
end interface

!------------------------------------------------------------------------------
