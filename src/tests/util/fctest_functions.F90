! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the logging facilities
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_logging_fxt
use atlas_module
!use atlas_functions_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg

end module fcta_logging_fxt

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_logging,fcta_logging_fxt)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_atlas_Functions )
  type(atlas_Functions) :: functions
  real(c_double) :: val_r8
  real(c_float) :: val_r4
  val_r8 = functions%MDPI_sinusoid(1._c_double, 1._c_double)
  val_r4 = functions%MDPI_harmonic(1., 1.)
  val_r4 = functions%MDPI_vortex(1., 1.)
  val_r8 = functions%MDPI_gulfstream(1., 1.)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
