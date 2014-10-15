! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: ReducedGG_type

! Purpose :
! -------
!   *ReducedGG* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: ngptot => ReducedGG__ngptot
  procedure :: nlat => ReducedGG__nlat
  procedure :: nlon => ReducedGG__nlon
  procedure :: lats => ReducedGG__lats
END TYPE ReducedGG_type

interface new_GaussianGrid
  module procedure new_reduced_gaussian_grid
  module procedure new_regular_gaussian_grid
  module procedure new_custom_reduced_gaussian_grid
end interface

interface new_LatlonGrid
  module procedure new_regular_latlon_grid
end interface

!------------------------------------------------------------------------------
