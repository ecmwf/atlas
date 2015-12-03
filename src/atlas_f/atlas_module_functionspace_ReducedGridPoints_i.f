! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_NextFunctionSpace) :: atlas_functionspace_ReducedGridPoints

! Purpose :
! -------
!   *atlas_functionspace_ReducedGridPoints* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: create_field_name     => ReducedGridPoints__create_field_name
  procedure, private :: create_field_name_lev => ReducedGridPoints__create_field_name_lev
  generic, public :: create_field => &
    & create_field_name, &
    & create_field_name_lev

  procedure, private :: create_glb_field_name     => ReducedGridPoints__create_glb_field_name
  procedure, private :: create_glb_field_name_lev => ReducedGridPoints__create_glb_field_name_lev
  generic, public :: create_global_field => &
    & create_glb_field_name, &
    & create_glb_field_name_lev

  procedure, public :: gather => ReducedGridPoints__gather
  procedure, public :: scatter => ReducedGridPoints__scatter
  
  procedure, private :: checksum_fieldset => ReducedGridPoints__checksum_fieldset
  procedure, private :: checksum_field => ReducedGridPoints__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field
  

#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_ReducedGridPoints__final
#endif

END TYPE atlas_functionspace_ReducedGridPoints

interface atlas_functionspace_ReducedGridPoints
  module procedure ReducedGridPoints__cptr
  module procedure ReducedGridPoints__grid
end interface


!------------------------------------------------------------------------------

