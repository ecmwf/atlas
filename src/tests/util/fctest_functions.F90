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
  FCTEST_CHECK_CLOSE(val, 1.0002115216773033_c_double, 1e-12_c_double)
  val = MDPI_harmonic(1._c_double, 1._c_double)
  FCTEST_CHECK_CLOSE(val, 2.0000000000000000_c_double, 1e-12_c_double)
  val = MDPI_vortex(1._c_double, 1._c_double)
  FCTEST_CHECK_CLOSE(val, 2.7267489215500755_c_double, 1e-12_c_double)
  val = MDPI_gulfstream(1._c_double, 1._c_double)
  FCTEST_CHECK_CLOSE(val, 1.0002115216773033_c_double, 1e-12_c_double)
END_TEST

TEST( test_atlas_Functions_vector )
  real(c_double), dimension(3) :: val, lon, lat, val_ref
  lon = [ 1._c_double, 2._c_double, 3._c_double ]
  lat = [ 1._c_double, 2._c_double, 3._c_double ]
  val = MDPI_sinusoid(lon, lat)
  val_ref = [1.0002115216773033_c_double, 1.0008458683590891_c_double, 1.0019023851484181_c_double]
  FCTEST_CHECK_CLOSE(val, val_ref, 1e-12_c_double)
  val = MDPI_harmonic(lon, lat)
  val_ref = [2.0000000000000000_c_double, 2.0000000000000000_c_double, 2.0000000000000000_c_double]
  FCTEST_CHECK_CLOSE(val, val_ref, 1e-12_c_double)
  val = MDPI_vortex(lon, lat)
  val_ref = [2.7267489215500755_c_double, 2.7520839004022091_c_double, 2.7755506683886928_c_double]
  FCTEST_CHECK_CLOSE(val, val_ref, 1e-12_c_double)
  val = MDPI_gulfstream(lon, lat)
  val_ref = [1.0002115216773033_c_double, 1.0008458683590891_c_double, 1.0019023851484181_c_double]
  FCTEST_CHECK_CLOSE(val, val_ref, 1e-12_c_double)
END_TEST

TEST( test_initialise_field )
  type(atlas_Field) :: field_xy, field_val
  real(c_double), dimension(:,:), pointer :: field_xy_v
  real(c_double), dimension(:), pointer :: field_val_v
  real(c_double), dimension(3) :: field_val_ref
  field_xy = atlas_Field(kind=atlas_real(c_double), shape=[2,3])
  field_val = atlas_Field(kind=atlas_real(c_double), shape=[3])
  call field_xy%data(field_xy_v)
  field_xy_v = 1._c_double
  call field_val%data(field_val_v)
  field_val_v = MDPI_sinusoid(field_xy_v(1,:), field_xy_v(2,:))
  field_val_ref = [1.0002115216773033_c_double, 1.0002115216773033_c_double, 1.0002115216773033_c_double]
  FCTEST_CHECK_CLOSE(field_val_v, field_val_ref, 1e-12_c_double)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
