! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module projection_fixture
use atlas_module
end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_projection,projection_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_rotated_schmidt_config )
type(atlas_Projection) :: projection
type(atlas_Config) :: config
config = atlas_Config()
call config%set("type","rotated_schmidt")
call config%set("stretching_factor",2.0_dp)
call config%set("rotation_angle", 180.0_dp)
call config%set("north_pole", [2.0_dp,46.7_dp] )
projection = atlas_Projection(config)
FCTEST_CHECK_EQUAL( projection%type(), "rotated_schmidt" )
FCTEST_CHECK_EQUAL( projection%hash(), "148d7ceb58250c0f48dc6b590941a341" )
call config%final()
call projection%final()
END_TEST

TEST( test_rotated_schmidt_specific_constructor )
type(atlas_Projection) :: projection
projection = atlas_RotatedSchmidtProjection( stretching_factor=2.0_dp, &
        north_pole=[2.0_dp,46.7_dp], rotation_angle=180.0_dp)
FCTEST_CHECK_EQUAL( projection%type(), "rotated_schmidt" )
FCTEST_CHECK_EQUAL( projection%hash(), "148d7ceb58250c0f48dc6b590941a341" )
call projection%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_lambert_conformal_conic_config )
type(atlas_Projection) :: projection
type(atlas_Config) :: config
config = atlas_Config()
call config%set("type","lambert_conformal_conic")
call config%set("longitude0",4.0_dp)
call config%set("latitude0", 50.0_dp)
projection = atlas_Projection(config)
FCTEST_CHECK_EQUAL( projection%type(), "lambert_conformal_conic" )
FCTEST_CHECK_EQUAL( projection%hash(), "47244ee950cf8357763a9b5b871b52ab" )
call config%final()
call projection%final()
END_TEST

TEST( test_lambert_conformal_conic_specific_constructor )
type(atlas_Projection) :: projection
projection = atlas_LambertConformalConicProjection(4.0_dp,50._dp)
FCTEST_CHECK_EQUAL( projection%type(), "lambert_conformal_conic" )
FCTEST_CHECK_EQUAL( projection%hash(), "47244ee950cf8357763a9b5b871b52ab" )
call projection%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_rotated_lonlat_config )
type(atlas_Projection) :: projection
type(atlas_Config) :: config
config = atlas_Config()
call config%set("type","rotated_lonlat")
call config%set("north_pole", [2.0_dp,46.7_dp] )
call config%set("rotation_angle", 180.0_dp)
projection = atlas_Projection(config)
FCTEST_CHECK_EQUAL( projection%type(), "rotated_lonlat" )
FCTEST_CHECK_EQUAL( projection%hash(), "2b6db0e1ccbe7c514dd726f408f92adb" )
call config%final()
call projection%final()
END_TEST

TEST( test_rotated_lonlat_specific_constructor )
type(atlas_Projection) :: projection
projection = atlas_RotatedLonLatProjection([2.0_dp,46.7_dp],180._dp)
FCTEST_CHECK_EQUAL( projection%type(), "rotated_lonlat" )
FCTEST_CHECK_EQUAL( projection%hash(), "2b6db0e1ccbe7c514dd726f408f92adb" )
call projection%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

