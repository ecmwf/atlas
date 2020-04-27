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

module fcta_UnitSphere_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_UnitSphere_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_UnitSphere,fcta_UnitSphere_fixture)

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

TEST( test_UnitSphere )
use fckit_log_module
use fckit_c_interop_module
implicit none

  type(atlas_UnitSphere) :: UnitSphere
  type(atlas_PointLonLat) :: pointLonLatA, pointLonLatB, pointLonLatA_test
  type(atlas_PointXYZ) :: pointXYZA, pointXYZB, pointXYZA_test
  real(c_double) :: Clat, Clon1, Clon2

  write(*,*) "test_UnitSphere starting"

  ! Check constructor
  UnitSphere = atlas_UnitSphere()
  write(0,*) "UnitSphere%c_ptr() = ", c_ptr_to_loc(UnitSphere%CPTR_PGIBUG_A)

  ! Check radius
  FCTEST_CHECK_EQUAL( UnitSphere%radius() , 1.0_c_double )

  ! Define points
  pointLonLatA = atlas_PointLonLat( -7.16e1_c_double, -3.3e1_c_double)
  pointLonLatB = atlas_PointLonLat( 1.218e2_c_double, 3.14e1_c_double)
  pointXYZA = atlas_PointXYZ( 2.647e-1_c_double, -7.958e-1_c_double, -5.446e-1_c_double)
  pointXYZB = atlas_PointXYZ( -4.498e-1_c_double, 7.254e-1_c_double, 5.210e-1_c_double)

  ! Check central angle
  FCTEST_CHECK_CLOSE( UnitSphere%central_angle(pointLonLatA, pointLonLatB) , 2.942_c_double , 1.e-3_c_double )
  FCTEST_CHECK_CLOSE( UnitSphere%central_angle(pointXYZA, pointXYZB) , 2.942_c_double , 1.e-3_c_double )

  ! Check distance
  FCTEST_CHECK_CLOSE( UnitSphere%distance(pointLonLatA, pointLonLatB) , 2.942_c_double , 1.e-3_c_double )
  FCTEST_CHECK_CLOSE( UnitSphere%distance(pointXYZA, pointXYZB) , 2.942_c_double , 1.e-3_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( UnitSphere%area() , 1.257e1_c_double , 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( UnitSphere%area(pointLonLatB, pointLonLatA), 3.099_c_double , 1.e-3_c_double )

  ! Check great circle
  Clat = UnitSphere%great_circle_latitude_given_longitude(pointLonLatA, pointLonLatB, -159.18_c_double)
  FCTEST_CHECK_CLOSE( Clat , -6.806_c_double , 1.e-3_c_double )
  call UnitSphere%great_circle_longitude_given_latitude(pointLonLatA, pointLonLatB, 0.0_c_double, Clon1, Clon2)
  FCTEST_CHECK_CLOSE( Clon1 , 370.3_c_double , 1.e-1_c_double )
  FCTEST_CHECK_CLOSE( Clon2 , 190.3_c_double , 1.e-1_c_double )

  ! Check conversion from spherical to cartesian
  pointXYZA_test = atlas_PointXYZ()
  call UnitSphere%convert_spherical_to_cartesian(pointLonLatA, pointXYZA_test)
  FCTEST_CHECK_CLOSE( pointXYZA_test%x() , pointXYZA%x(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( pointXYZA_test%y() , pointXYZA%y(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( pointXYZA_test%z() , pointXYZA%z(), 1.e3_c_double )

  ! Check conversion from cartesian to spherical
  pointLonLatA_test = atlas_PointLonLat()
  call UnitSphere%convert_cartesian_to_spherical(pointXYZA, pointLonLatA_test)
  FCTEST_CHECK_CLOSE( pointLonLatA_test%lon() , pointLonLatA%lon(), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( pointLonLatA_test%lat() , pointLonLatA%lat(), 1.e-2_c_double )

  call UnitSphere%final()
  call pointLonLatA%final()
  call pointLonLatB%final()
  call pointXYZA%final()
  call pointXYZB%final()
END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

