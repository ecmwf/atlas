! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_NextFunctionSpace) :: atlas_functionspace_ReducedGridPoint

! Purpose :
! -------
!   *atlas_functionspace_ReducedGridPoint* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: create_field_name     => ReducedGridPoint__create_field_name
  procedure, private :: create_field_name_lev => ReducedGridPoint__create_field_name_lev
  generic, public :: create_field => &
    & create_field_name, &
    & create_field_name_lev

  procedure, private :: create_glb_field_name     => ReducedGridPoint__create_glb_field_name
  procedure, private :: create_glb_field_name_lev => ReducedGridPoint__create_glb_field_name_lev
  generic, public :: create_global_field => &
    & create_glb_field_name, &
    & create_glb_field_name_lev

  procedure, public :: gather => ReducedGridPoint__gather
  procedure, public :: scatter => ReducedGridPoint__scatter
  
  procedure, private :: checksum_fieldset => ReducedGridPoint__checksum_fieldset
  procedure, private :: checksum_field => ReducedGridPoint__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field
  

#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_ReducedGridPoint__final
#endif

END TYPE atlas_functionspace_ReducedGridPoint

interface atlas_functionspace_ReducedGridPoint
  module procedure ReducedGridPoint__cptr
  module procedure ReducedGridPoint__grid
end interface


!------------------------------------------------------------------------------

