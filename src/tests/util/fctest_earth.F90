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

module fcta_Earth_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_Earth_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Earth,fcta_Earth_fixture)

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

TEST( test_earth )
use fckit_log_module
use fckit_c_interop_module
implicit none

  type(atlas_Earth) :: earth
  type(atlas_PointLonLat) :: pointLonLatA, pointLonLatB, pointLonLatA_test
  type(atlas_PointXYZ) :: pointXYZA, pointXYZB, pointXYZA_test
  real(c_double) :: Clat, Clon1, Clon2

  write(*,*) "test_earth starting"

  ! Check constructor
  earth = atlas_Earth()
  write(0,*) "earth%c_ptr() = ", c_ptr_to_loc(earth%CPTR_PGIBUG_A)

  ! Check radius
  FCTEST_CHECK_EQUAL( earth%radius() , 6371229.0_c_double )

  ! Define points
  pointLonLatA = atlas_PointLonLat( -71.6_c_double, -33.0_c_double)
  pointLonLatB = atlas_PointLonLat( 121.8_c_double, 31.4_c_double)
  pointXYZA = atlas_PointXYZ( 1.687e6_c_double, -5.070e6_c_double, -3.470e6_c_double)
  pointXYZB = atlas_PointXYZ( -2.866e6_c_double, 4.622e6_c_double, 3.319e6_c_double)

  ! Check central angle
  FCTEST_CHECK_CLOSE( earth%central_angle(pointLonLatA, pointLonLatB) , 2.942_c_double , 1.e-3_c_double )
  FCTEST_CHECK_CLOSE( earth%central_angle(pointXYZA, pointXYZB) , 2.942_c_double , 1.e-3_c_double )

  ! Check distance
  FCTEST_CHECK_CLOSE( earth%distance(pointLonLatA, pointLonLatB) , 1.874e7_c_double , 1.e4_c_double )
  FCTEST_CHECK_CLOSE( earth%distance(pointXYZA, pointXYZB) , 1.874e7_c_double , 1.e4_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( earth%area() , 5.101e14_c_double , 1.e11_c_double )
  FCTEST_CHECK_CLOSE( earth%area(pointLonLatB, pointLonLatA), 1.258e14_c_double , 1.e11_c_double )

  ! Check great circle
  Clat = earth%great_circle_latitude_given_longitude(pointLonLatA, pointLonLatB, -159.18_c_double)
  FCTEST_CHECK_CLOSE( Clat , -6.806_c_double , 1.e-3_c_double )
  call earth%great_circle_longitude_given_latitude(pointLonLatA, pointLonLatB, 0.0_c_double, Clon1, Clon2)
  FCTEST_CHECK_CLOSE( Clon1 , 370.3_c_double , 1.e-1_c_double )
  FCTEST_CHECK_CLOSE( Clon2 , 190.3_c_double , 1.e-1_c_double )

  ! Check conversion from spherical to cartesian
  pointXYZA_test = atlas_PointXYZ()
  call earth%convert_spherical_to_cartesian(pointLonLatA, pointXYZA_test)
  FCTEST_CHECK_CLOSE( pointXYZA_test%x() , pointXYZA%x(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( pointXYZA_test%y() , pointXYZA%y(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( pointXYZA_test%z() , pointXYZA%z(), 1.e3_c_double )

  ! Check conversion from cartesian to spherical
  pointLonLatA_test = atlas_PointLonLat()
  call earth%convert_cartesian_to_spherical(pointXYZA, pointLonLatA_test)
  FCTEST_CHECK_CLOSE( pointLonLatA_test%lon() , pointLonLatA%lon(), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( pointLonLatA_test%lat() , pointLonLatA%lat(), 1.e-2_c_double )

  call earth%final()
  call pointLonLatA%final()
  call pointLonLatB%final()
  call pointXYZA%final()
  call pointXYZB%final()
END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

