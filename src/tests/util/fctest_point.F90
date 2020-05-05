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

  type(atlas_PointXY) :: pointXY
  type(atlas_PointXYZ) :: pointXYZ
  type(atlas_PointLonLat) :: pointLonLat
  type(fckit_logchannel) :: info

  ! Get info channel
  info = fckit_log%info_channel()

  write(*,*) "test_point_xy starting"

  ! Check empty constructor and assign
  pointXY = atlas_PointXY()
  write(0,*) "pointXY%c_ptr() = ", c_ptr_to_loc(pointXY%CPTR_PGIBUG_A)
  call pointXY%assign(1.687e6_c_double, -5.070e6_c_double)
  FCTEST_CHECK_EQUAL( pointXY%x() , 1.687e6_c_double )
  FCTEST_CHECK_EQUAL( pointXY%y() , -5.070e6_c_double )
  call pointXY%final()

  ! Check constructor with x/y
  pointXY = atlas_PointXY(1.687e6_c_double, -5.070e6_c_double)
  FCTEST_CHECK_EQUAL( pointXY%x() , 1.687e6_c_double )
  FCTEST_CHECK_EQUAL( pointXY%y() , -5.070e6_c_double )

  ! Check print
  call pointXY%print(info)

  ! Finalization
  call pointXY%final()

  write(*,*) "test_point_xyz starting"

  ! Check empty constructor and assign
  pointXYZ = atlas_PointXYZ()
  write(0,*) "pointXYZ%c_ptr() = ", c_ptr_to_loc(pointXYZ%CPTR_PGIBUG_A)
  call pointXYZ%assign(1.687e6_c_double, -5.070e6_c_double, 3.470e6_c_double)
  FCTEST_CHECK_EQUAL( pointXYZ%x() , 1.687e6_c_double )
  FCTEST_CHECK_EQUAL( pointXYZ%y() , -5.070e6_c_double )
  FCTEST_CHECK_EQUAL( pointXYZ%z() , 3.470e6_c_double )
  call pointXYZ%final()

  ! Check constructor with x/y/z
  pointXYZ = atlas_PointXYZ(1.687e6_c_double, -5.070e6_c_double, 3.470e6_c_double)
  FCTEST_CHECK_EQUAL( pointXYZ%x() , 1.687e6_c_double )
  FCTEST_CHECK_EQUAL( pointXYZ%y() , -5.070e6_c_double )
  FCTEST_CHECK_EQUAL( pointXYZ%z() , 3.470e6_c_double )

  ! Check print
  call pointXYZ%print(info)

  ! Finalization
  call pointXYZ%final()

  write(*,*) "test_point_lonlat starting"

  ! Check empty constructor and assign
  pointLonLat = atlas_PointLonLat()
  write(0,*) "pointLonLat%c_ptr() = ", c_ptr_to_loc(pointLonLat%CPTR_PGIBUG_A)
  call pointLonLat%assign(-71.6_c_double, -33._c_double)
  FCTEST_CHECK_EQUAL( pointLonLat%lon() , -71.6_c_double )
  FCTEST_CHECK_EQUAL( pointLonLat%lat() , -33._c_double )
  call pointLonLat%final()

  ! Check constructor with lon/lat
  pointLonLat = atlas_PointLonLat(-71.6_c_double, -33._c_double)
  FCTEST_CHECK_EQUAL( pointLonLat%lon() , -71.6_c_double )
  FCTEST_CHECK_EQUAL( pointLonLat%lat() , -33._c_double )

  ! Check print
  call pointLonLat%print(info)

  ! Check normalise
  call pointLonLat%normalise()
  FCTEST_CHECK_CLOSE( pointLonLat%lon() , -71.6_c_double + 360.0_c_double, 1.e-12_c_double )
  call pointLonLat%normalise(-180.0_c_double, 180.0_c_double)
  FCTEST_CHECK_CLOSE( pointLonLat%lon() , -71.6_c_double, 1.e-12_c_double )

  ! Finalization
  call pointLonLat%final()

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

