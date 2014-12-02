! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: ReducedGrid_type

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
END TYPE ReducedGrid_type

interface new_GaussianGrid
  module procedure new_reduced_gaussian_grid
  module procedure new_regular_gaussian_grid
  module procedure new_custom_reduced_gaussian_grid
end interface

interface new_LatlonGrid
  module procedure new_regular_latlon_grid
end interface

!------------------------------------------------------------------------------
