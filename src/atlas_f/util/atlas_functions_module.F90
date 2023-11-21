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
use atlas_functions_c_binding

implicit none

contains

function MDPI_sinusoid(lon, lat) result(val)
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_sinusoid(lon, lat)
end function MDPI_sinusoid

! -----------------------------------------------------------------------------

function MDPI_harmonic(lon, lat) result(val)
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_harmonic(lon, lat)
end function MDPI_harmonic

! -----------------------------------------------------------------------------

function MDPI_vortex(lon, lat) result(val)
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_vortex(lon, lat)
end function MDPI_vortex

! -----------------------------------------------------------------------------

function MDPI_gulfstream(lon, lat) result(val)
  real(c_double), intent(in) :: lon, lat
  real(c_double) :: val
  val = atlas__functions__MDPI_gulfstream(lon, lat)
end function MDPI_gulfstream

! -----------------------------------------------------------------------------

end module atlas_functions_module
