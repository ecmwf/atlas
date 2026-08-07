! (C) Copyright 2026 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_algorithms_spectral_module

use, intrinsic :: iso_c_binding, only : c_int, c_double
use atlas_Field_module, only : atlas_Field
use atlas_functionspace_Spectral_module, only : atlas_functionspace_Spectral

implicit none

private :: c_int, c_double
private :: atlas_Field
private :: atlas_functionspace_Spectral

public :: filter_spectral_cutoff
public :: power_spectrum

interface power_spectrum
  module procedure power_spectrum_rank1
  module procedure power_spectrum_rank2
end interface

contains

subroutine filter_spectral_cutoff(spectral, field, cutoff)
  use atlas_algorithms_spectral_filter_cutoff_c_binding
  type(atlas_functionspace_Spectral), intent(in) :: spectral
  type(atlas_Field), intent(inout) :: field
  integer(c_int), intent(in) :: cutoff

  call atlas__spectral__filter_cutoff(spectral%c_ptr(), field%c_ptr(), cutoff)
end subroutine filter_spectral_cutoff

! -----------------------------------------------------------------------------

subroutine power_spectrum_rank1(spectral, field, spectrum)
  use atlas_algorithms_spectral_power_spectrum_c_binding
  type(atlas_functionspace_Spectral), intent(in) :: spectral
  type(atlas_Field), intent(in) :: field
  real(c_double), intent(inout) :: spectrum(:)
  integer(c_int) :: spectrum_extents(2)

  spectrum_extents = [int(size(spectrum), c_int), 1_c_int]
  call atlas__spectral__power_spectrum(spectral%c_ptr(), field%c_ptr(), spectrum, spectrum_extents)
end subroutine power_spectrum_rank1

! -----------------------------------------------------------------------------

subroutine power_spectrum_rank2(spectral, field, spectrum)
  use atlas_algorithms_spectral_power_spectrum_c_binding
  type(atlas_functionspace_Spectral), intent(in) :: spectral
  type(atlas_Field), intent(in) :: field
  real(c_double), intent(inout) :: spectrum(:,:)
  integer(c_int) :: spectrum_extents(2)

  spectrum_extents = [int(size(spectrum, 2), c_int), int(size(spectrum, 1), c_int)]
  call atlas__spectral__power_spectrum(spectral%c_ptr(), field%c_ptr(), spectrum, spectrum_extents)
end subroutine power_spectrum_rank2

! -----------------------------------------------------------------------------

end module atlas_algorithms_spectral_module
