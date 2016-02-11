! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_Grid

! Purpose :
! -------
!   *atlas_Grid* : Object Grid specifications for Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: npts => atlas_Grid__npts
  procedure, public :: delete => atlas_Grid__delete
  procedure, public :: copy => atlas_Grid__copy
END TYPE atlas_Grid

!------------------------------------------------------------------------------

TYPE, extends(atlas_Grid) :: atlas_ReducedGrid

! Purpose :
! -------
!   *atlas_ReducedGrid* : Object Grid specifications for Reduced Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: N    => ReducedGrid__N
  procedure :: nlat => ReducedGrid__nlat
  procedure :: nlon => ReducedGrid__nlon
  procedure :: nlonmax => ReducedGrid__nlonmax
  procedure :: lat => ReducedGrid__latitudes
END TYPE atlas_ReducedGrid

!------------------------------------------------------------------------------

TYPE, extends(atlas_ReducedGrid) :: atlas_ReducedGaussianGrid

! Purpose :
! -------
!   *atlas_ReducedGaussianGrid* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_ReducedGaussianGrid

!------------------------------------------------------------------------------

TYPE, extends(atlas_ReducedGrid) :: atlas_GaussianGrid

! Purpose :
! -------
!   *atlas_GaussianGrid* : Object Grid specifications for Regular Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_GaussianGrid

!------------------------------------------------------------------------------

TYPE, extends(atlas_ReducedGrid) :: atlas_LonLatGrid

! Purpose :
! -------
!   *atlas_LonLatGrid* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_LonLatGrid

!------------------------------------------------------------------------------

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
