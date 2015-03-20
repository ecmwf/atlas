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
  procedure :: handle       => atlas_Trans__handle
  procedure :: nproc        => atlas_Trans__nproc
  procedure :: myproc       => atlas_Trans__myproc
  procedure :: ndgl         => atlas_Trans__ndgl
  procedure :: nsmax        => atlas_Trans__nsmax
  procedure :: ngptot       => atlas_Trans__ngptot
  procedure :: ngptotg      => atlas_Trans__ngptotg
  procedure :: ngptotmx     => atlas_Trans__ngptotmx
  procedure :: nspec        => atlas_Trans__nspec
  procedure :: nspec2       => atlas_Trans__nspec2
  procedure :: nspec2g      => atlas_Trans__nspec2g
  procedure :: nspec2mx     => atlas_Trans__nspec2mx
  procedure :: n_regions_NS => atlas_Trans__n_regions_NS
  procedure :: n_regions_EW => atlas_Trans__n_regions_EW
  procedure :: nloen        => atlas_Trans__nloen
  procedure :: n_regions    => atlas_Trans__n_regions
  procedure :: nfrstlat     => atlas_Trans__nfrstlat
  procedure :: nlstlat      => atlas_Trans__nlstlat
  procedure :: nptrfrstlat  => atlas_Trans__nptrfrstlat
  procedure :: nsta         => atlas_Trans__nsta
  procedure :: nonl         => atlas_Trans__nonl
  procedure :: nmyms        => atlas_Trans__nmyms
  procedure :: nasm0        => atlas_Trans__nasm0
  procedure :: nump         => atlas_Trans__nump

END TYPE atlas_Trans

!------------------------------------------------------------------------------

interface new_atlas_Trans
  module procedure new_atlas_Trans
end interface
