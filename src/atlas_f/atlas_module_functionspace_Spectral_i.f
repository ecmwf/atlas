! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_Spectral

! Purpose :
! -------
!   *atlas_functionspace_Spectral* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: create_field_name     => SpectralFunctionSpace__create_field_name
  procedure, private :: create_field_name_lev => SpectralFunctionSpace__create_field_name_lev
  generic, public :: create_field => &
    & create_field_name, &
    & create_field_name_lev

  procedure, private :: create_glb_field_name     => SpectralFunctionSpace__create_glb_field_name
  procedure, private :: create_glb_field_name_lev => SpectralFunctionSpace__create_glb_field_name_lev
  generic, public :: create_global_field => &
    & create_glb_field_name, &
    & create_glb_field_name_lev

  procedure, public :: gather => SpectralFunctionspace__gather
  procedure, public :: scatter => SpectralFunctionspace__scatter

#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_Spectral__final
#endif

END TYPE atlas_functionspace_Spectral

interface atlas_functionspace_Spectral
  module procedure atlas_functionspace_Spectral__cptr
  module procedure atlas_functionspace_Spectral__truncation
  module procedure atlas_functionspace_Spectral__trans
end interface


!------------------------------------------------------------------------------

