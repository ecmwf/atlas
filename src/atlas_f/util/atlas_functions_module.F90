! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_functions_module

use, intrinsic :: iso_c_binding
use fckit_owned_object_module, only : fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_functions

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_functions

! Purpose :
! -------
!   *function* : Analytic functions for idealised initial data

! Methods :
! -------

! Author :
! ------
!   Nov 2023 Slavko Brdar (ECMWF)

!------------------------------------------------------------------------------
contains
  procedure, public :: MDPI_sinusoid_r8 => atlas_functions__MDPI_sinusoid_r8
  procedure, public :: MDPI_sinusoid_r4 => atlas_functions__MDPI_sinusoid_r4
  procedure, public :: MDPI_harmonic_r8 => atlas_functions__MDPI_harmonic_r8
  procedure, public :: MDPI_harmonic_r4 => atlas_functions__MDPI_harmonic_r4
  procedure, public :: MDPI_vortex_r8 => atlas_functions__MDPI_vortex_r8
  procedure, public :: MDPI_vortex_r4 => atlas_functions__MDPI_vortex_r4
  procedure, public :: MDPI_gulfstream_r8 => atlas_functions__MDPI_gulfstream_r8
  procedure, public :: MDPI_gulfstream_r4 => atlas_functions__MDPI_gulfstream_r4
  generic :: MDPI_sinusoid => MDPI_sinusoid_r8, MDPI_sinusoid_r4
  generic :: MDPI_harmonic => MDPI_harmonic_r8, MDPI_harmonic_r4
  generic :: MDPI_vortex => MDPI_vortex_r8, MDPI_vortex_r4
  generic :: MDPI_gulfstream => MDPI_gulfstream_r8, MDPI_gulfstream_r4
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_functions__final_auto
#endif

END TYPE atlas_functions

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------

function atlas_functions__MDPI_sinusoid_r8(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_sinusoid_r8(lon, lat)
end function atlas_functions__MDPI_sinusoid_r8

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_harmonic_r8(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_harmonic_r8(lon, lat)
end function atlas_functions__MDPI_harmonic_r8

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_vortex_r8(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_vortex_r8(lon, lat)
end function atlas_functions__MDPI_vortex_r8

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_gulfstream_r8(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_gulfstream_r8(lon, lat)
end function atlas_functions__MDPI_gulfstream_r8

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_sinusoid_r4(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_float), intent(in) :: lon, lat
  real(c_float) :: val
  val = atlas__functions__MDPI_sinusoid_r4(lon, lat)
end function atlas_functions__MDPI_sinusoid_r4

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_harmonic_r4(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_float), intent(in) :: lon, lat
  real(c_float) :: val
  val = atlas__functions__MDPI_harmonic_r4(lon, lat)
end function atlas_functions__MDPI_harmonic_r4

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_vortex_r4(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_float), intent(in) :: lon, lat
  real(c_float) :: val
  val = atlas__functions__MDPI_vortex_r4(lon, lat)
end function atlas_functions__MDPI_vortex_r4

! -----------------------------------------------------------------------------

function atlas_functions__MDPI_gulfstream_r4(this, lon, lat) result(val)
  use atlas_functions_c_binding
  class(atlas_functions), intent(in) :: this
  real(c_float), intent(in) :: lon, lat
  real(c_float) :: val
  val = atlas__functions__MDPI_gulfstream_r4(lon, lat)
end function atlas_functions__MDPI_gulfstream_r4

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_functions__final_auto(this)
  type(atlas_functions), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_function__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_functions_module
