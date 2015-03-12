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

interface new_atlas_ReducedGrid
  module procedure new_atlas_reduced_grid
  module procedure new_atlas_gaussian_grid
  module procedure new_atlas_reduced_gaussian_grid
  module procedure new_atlas_lonlat_grid
end interface

!------------------------------------------------------------------------------
