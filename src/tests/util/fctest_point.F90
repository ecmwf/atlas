! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"
#include "atlas/atlas_f.h"

! -----------------------------------------------------------------------------

module fcta_Point_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_Point_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Point,fcta_Point_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use fckit_main_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_point )
use fckit_log_module
use fckit_c_interop_module
implicit none

  type(atlas_PointLonLat) :: pointLonLat
  type(fckit_logchannel) :: info

  write(*,*) "test_point starting"

  ! Check empty constructor and assign
  pointLonLat = atlas_PointLonLat()
  write(0,*) "pointLonLat%c_ptr() = ", c_ptr_to_loc(pointLonLat%CPTR_PGIBUG_A)
  call pointLonLat%assign(-71.6_c_double, -33.0_c_double)
  FCTEST_CHECK_EQUAL( pointLonLat%lon() , -71.6_c_double )
  FCTEST_CHECK_EQUAL( pointLonLat%lat() , -33.0_c_double )
  call pointLonLat%final()

  ! Check constructor with lon/lat
  pointLonLat = atlas_PointLonLat(-71.6_c_double, -33.0_c_double)
  FCTEST_CHECK_EQUAL( pointLonLat%lon() , -71.6_c_double )
  FCTEST_CHECK_EQUAL( pointLonLat%lat() , -33.0_c_double )

  ! Check print
  info = fckit_log%info_channel()
  call pointLonLat%print(info)

  ! Check normalise
  call pointLonLat%normalise()
  FCTEST_CHECK_CLOSE( pointLonLat%lon() , -71.6_c_double + 360.0_c_double, 1.e-12_c_double )
  call pointLonLat%normalise(-180.0_c_double, 180.0_c_double)
  FCTEST_CHECK_CLOSE( pointLonLat%lon() , -71.6_c_double, 1.e-12_c_double )
  call pointLonLat%final()

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

