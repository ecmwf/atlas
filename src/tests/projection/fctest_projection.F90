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

TEST( test_xy2lonlat )
type(atlas_Projection) :: projection
type(atlas_Config) :: config
real(dp) :: x, y, lon, lat
config = atlas_Config()
call config%set("type","mercator")
call config%set("latitude1",45.0_dp)
projection = atlas_Projection(config)
x = 5000.0e3_dp
y = 5000.0e3_dp
call projection%xy2lonlat(x, y, lon, lat)
FCTEST_CHECK_CLOSE( lon, 63.589_dp, 1.0e-3_dp )
FCTEST_CHECK_CLOSE( lat, 53.514_dp, 1.0e-3_dp )
call config%final()
call projection%final()
END_TEST

TEST( test_lonlat2xy )
type(atlas_Projection) :: projection
type(atlas_Config) :: config
real(dp) :: lon, lat, x, y
config = atlas_Config()
call config%set("type","mercator")
call config%set("latitude1",45.0_dp)
projection = atlas_Projection(config)
lon = 12.0_dp
lat = 38.0_dp
call projection%lonlat2xy(lon, lat, x, y)
FCTEST_CHECK_CLOSE( x, 943554.154_dp, 1.0e-3_dp )
FCTEST_CHECK_CLOSE( y, 3234635.895_dp, 1.0e-3_dp )
call config%final()
call projection%final()
END_TEST

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
FCTEST_CHECK_EQUAL( projection%hash(), "3a39e0635b7d0a45f684696ca89825e6" )
call config%final()
call projection%final()
END_TEST

TEST( test_rotated_schmidt_specific_constructor )
type(atlas_Projection) :: projection
projection = atlas_RotatedSchmidtProjection( stretching_factor=2.0_dp, &
        north_pole=[2.0_dp,46.7_dp], rotation_angle=180.0_dp)
FCTEST_CHECK_EQUAL( projection%type(), "rotated_schmidt" )
FCTEST_CHECK_EQUAL( projection%hash(), "3a39e0635b7d0a45f684696ca89825e6" )
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
FCTEST_CHECK_EQUAL( projection%hash(), "79586cfbc8145cdef1a25d075a9ae07e" )
call config%final()
call projection%final()
END_TEST

TEST( test_rotated_lonlat_specific_constructor )
type(atlas_Projection) :: projection
projection = atlas_RotatedLonLatProjection([2.0_dp,46.7_dp],180._dp)
FCTEST_CHECK_EQUAL( projection%type(), "rotated_lonlat" )
FCTEST_CHECK_EQUAL( projection%hash(), "79586cfbc8145cdef1a25d075a9ae07e" )
call projection%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

