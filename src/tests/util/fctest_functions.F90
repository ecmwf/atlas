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

module fcta_functions_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg

end module fcta_functions_fxt

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_functions,fcta_functions_fxt)

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
  real(c_double) :: val
  val = MDPI_sinusoid(1._c_double, 1._c_double)
  val = MDPI_harmonic(1._c_double, 1._c_double)
  val = MDPI_vortex(1._c_double, 1._c_double)
  val = MDPI_gulfstream(1._c_double, 1._c_double)
END_TEST

TEST( test_atlas_Functions_vector )
  real(c_double), dimension(3) :: val, lon, lat
  lon = [ 1._c_double, 1._c_double, 1._c_double ]
  lat = [ 1._c_double, 1._c_double, 1._c_double ]
  val = MDPI_sinusoid(lon, lat)
  val = MDPI_harmonic(lon, lat)
  val = MDPI_vortex(lon, lat)
  val = MDPI_gulfstream(lon, lat)
END_TEST

TEST( test_initialise_field )
  type(atlas_Field) :: field
  real(c_double), dimension(:,:), pointer :: fieldv
  real(c_double), dimension(3) :: val
  field = atlas_Field(kind=atlas_real(c_double), shape=[3,3])
  call field%data(fieldv)
  val = MDPI_sinusoid(fieldv(1,:), fieldv(2,:))
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
