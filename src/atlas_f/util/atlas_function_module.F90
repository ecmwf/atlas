! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_function_module

use, intrinsic :: iso_c_binding
use fckit_owned_object_module, only : fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_function

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_function

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
  procedure, public :: MDPI_sinusoid => atlas_function__MDPI_sinusoid
  procedure, public :: MDPI_harmonic => atlas_function__MDPI_harmonic
  procedure, public :: MDPI_vortex => atlas_function__MDPI_vortex
  procedure, public :: MDPI_gulfstream => atlas_function__MDPI_gulfstream
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_function__final_auto
#endif

END TYPE atlas_function

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------

function atlas_function__MDPI_sinusoid(this, lon, lat) result(val)
  use atlas_function_c_binding
  class(atlas_function), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__function__MDPI_sinusoid(lon, lat)
end function atlas_function__MDPI_sinusoid

! -----------------------------------------------------------------------------

function atlas_function__MDPI_harmonic(this, lon, lat) result(val)
  use atlas_function_c_binding
  class(atlas_function), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__function__MDPI_harmonic(lon, lat)
end function atlas_function__MDPI_harmonic

! -----------------------------------------------------------------------------

function atlas_function__MDPI_vortex(this, lon, lat) result(val)
  use atlas_function_c_binding
  class(atlas_function), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__function__MDPI_vortex(lon, lat)
end function atlas_function__MDPI_vortex

! -----------------------------------------------------------------------------

function atlas_function__MDPI_gulfstream(this, lon, lat) result(val)
  use atlas_function_c_binding
  class(atlas_function), intent(in) :: this
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__function__MDPI_gulfstream(lon, lat)
end function atlas_function__MDPI_gulfstream

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_function__final_auto(this)
  type(atlas_function), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_function__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_function_module
