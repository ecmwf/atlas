! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_SpectralFunctionSpace

! Purpose :
! -------
!   *atlas_SpectralFunctionSpace* : Interpretes spectral fields

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

  procedure :: finalize => atlas_SpectralFunctionSpace__finalize
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_SpectralFunctionSpace__final
#endif

END TYPE atlas_SpectralFunctionSpace

interface atlas_SpectralFunctionSpace
  module procedure atlas_SpectralFunctionSpace__cptr
  module procedure atlas_SpectralFunctionSpace__name_truncation
  module procedure atlas_SpectralFunctionSpace__name_trans
end interface


interface assignment(=)
  module procedure atlas_SpectralFunctionSpace__reset
end interface

!------------------------------------------------------------------------------

