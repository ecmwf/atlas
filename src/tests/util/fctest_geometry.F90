! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Geometry Datastructure
! @author Benjamin Menetrier

#include "fckit/fctest.h"
#include "atlas/atlas_f.h"

! -----------------------------------------------------------------------------

module fcta_Geometry_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_Geometry_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Geometry,fcta_Geometry_fixture)

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

TEST( test_geometry )
use fckit_log_module
use fckit_c_interop_module
implicit none

  type(atlas_Geometry) :: geometry
  real(c_double) :: lonlat1(2), lonlat2(2), lonlat(2)
  real(c_double) :: xyz1(3), xyz2(3), xyz(3)
  real(c_double) :: distance

  write(*,*) "test_geometry for UnitSphere starting"

  ! Check constructor
  geometry = atlas_Geometry("UnitSphere")
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Define points
  lonlat1 = (/-71.6_c_double, -33._c_double/)
  lonlat2 = (/121.8_c_double, 31.4_c_double/)
  xyz1 = (/2.647e-1_c_double, -7.958e-1_c_double, -5.446e-1_c_double/)
  xyz2 = (/-4.498e-1_c_double, 7.254e-1_c_double, 5.210e-1_c_double/)

  ! Check conversion from cartesian to spherical
  call geometry%xyz2lonlat(xyz1(1), xyz1(2), xyz1(3), lonlat(1), lonlat(2))
  FCTEST_CHECK_CLOSE( lonlat(1) , lonlat1(1), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lonlat(2) , lonlat1(2), 1.e-2_c_double )
  call geometry%xyz2lonlat(xyz1, lonlat)
  FCTEST_CHECK_CLOSE( lonlat(1) , lonlat1(1), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lonlat(2) , lonlat1(2), 1.e-2_c_double )

  ! Check conversion from spherical to cartesian
  call geometry%lonlat2xyz(lonlat1(1), lonlat1(2), xyz(1), xyz(2), xyz(3))
  FCTEST_CHECK_CLOSE( xyz(1) , xyz1(1) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(2) , xyz1(2) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(3) , xyz1(3) , 1.e3_c_double )
  call geometry%lonlat2xyz(lonlat1, xyz)
  FCTEST_CHECK_CLOSE( xyz(1) , xyz1(1) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(2) , xyz1(2) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(3) , xyz1(3) , 1.e3_c_double )

  ! Check distance
  distance = geometry%distance(lonlat1(1), lonlat1(2), lonlat2(1), lonlat2(2))
  FCTEST_CHECK_CLOSE( distance , 2.942_c_double , 1.e-3_c_double )
  distance = geometry%distance(lonlat1, lonlat2)
  FCTEST_CHECK_CLOSE( distance , 2.942_c_double , 1.e-3_c_double )
  distance = geometry%distance(xyz1(1), xyz1(2), xyz1(3), xyz2(1), xyz2(2), xyz2(3))
  FCTEST_CHECK_CLOSE( distance , 2.942_c_double , 1.e-3_c_double )
  distance = geometry%distance(xyz1, xyz2)
  FCTEST_CHECK_CLOSE( distance , 2.942_c_double , 1.e-3_c_double )

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry%radius() , 1.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry%area() , 1.257e1_c_double , 1.e-2_c_double )

  ! Finalization
  call geometry%final()

  write(*,*) "test_geometry for Earth starting"

  ! Check constructor for Earth
  geometry = atlas_Geometry("Earth")
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Define points
  lonlat1 = (/-71.6_c_double, -33._c_double/)
  lonlat2 = (/121.8_c_double, 31.4_c_double/)
  xyz1 = (/1.687e6_c_double, -5.070e6_c_double, -3.470e6_c_double/)
  xyz2 = (/-2.866e6_c_double, 4.622e6_c_double, 3.319e6_c_double/)

  ! Check conversion from cartesian to spherical
  call geometry%xyz2lonlat(xyz1(1), xyz1(2), xyz1(3), lonlat(1), lonlat(2))
  FCTEST_CHECK_CLOSE( lonlat(1) , lonlat1(1), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lonlat(2) , lonlat1(2), 1.e-2_c_double )
  call geometry%xyz2lonlat(xyz1, lonlat)
  FCTEST_CHECK_CLOSE( lonlat(1) , lonlat1(1), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lonlat(2) , lonlat1(2), 1.e-2_c_double )

  ! Check conversion from spherical to cartesian
  call geometry%lonlat2xyz(lonlat1(1), lonlat1(2), xyz(1), xyz(2), xyz(3))
  FCTEST_CHECK_CLOSE( xyz(1) , xyz1(1) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(2) , xyz1(2) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(3) , xyz1(3) , 1.e3_c_double )
  call geometry%lonlat2xyz(lonlat1, xyz)
  FCTEST_CHECK_CLOSE( xyz(1) , xyz1(1) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(2) , xyz1(2) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(3) , xyz1(3) , 1.e3_c_double )

  ! Check distance
  distance = geometry%distance(lonlat1(1), lonlat1(2), lonlat2(1), lonlat2(2))
  FCTEST_CHECK_CLOSE( distance , 1.874e7_c_double , 1.e4_c_double )
  distance = geometry%distance(lonlat1, lonlat2)
  FCTEST_CHECK_CLOSE( distance , 1.874e7_c_double , 1.e4_c_double )
  distance = geometry%distance(xyz1(1), xyz1(2), xyz1(3), xyz2(1), xyz2(2), xyz2(3))
  FCTEST_CHECK_CLOSE( distance , 1.874e7_c_double , 1.e4_c_double )
  distance = geometry%distance(xyz1, xyz2)
  FCTEST_CHECK_CLOSE( distance , 1.874e7_c_double , 1.e4_c_double )

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry%radius() , 6371229.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry%area() , 5.101e14_c_double , 1.e11_c_double )

  ! Finalization
  call geometry%final()

  write(*,*) "test_geometry for another planet (Mars here) starting"

  ! Check constructor for another planet (Mars here)
  geometry = atlas_Geometry( 3389500.0_c_double )
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Define points
  lonlat1 = (/-71.6_c_double, -33._c_double/)
  lonlat2 = (/121.8_c_double, 31.4_c_double/)
  xyz1 = (/8.973e5_c_double, -2.697e6_c_double, -1.846e6_c_double/)
  xyz2 = (/-1.525e6_c_double, 2.459e6_c_double, 1.766e6_c_double/)

  ! Check conversion from cartesian to spherical
  call geometry%xyz2lonlat(xyz1(1), xyz1(2), xyz1(3), lonlat(1), lonlat(2))
  FCTEST_CHECK_CLOSE( lonlat(1) , lonlat1(1), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lonlat(2) , lonlat1(2), 1.e-2_c_double )
  call geometry%xyz2lonlat(xyz1, lonlat)
  FCTEST_CHECK_CLOSE( lonlat(1) , lonlat1(1), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lonlat(2) , lonlat1(2), 1.e-2_c_double )

  ! Check conversion from spherical to cartesian
  call geometry%lonlat2xyz(lonlat1(1), lonlat1(2), xyz(1), xyz(2), xyz(3))
  FCTEST_CHECK_CLOSE( xyz(1) , xyz1(1) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(2) , xyz1(2) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(3) , xyz1(3) , 1.e3_c_double )
  call geometry%lonlat2xyz(lonlat1, xyz)
  FCTEST_CHECK_CLOSE( xyz(1) , xyz1(1) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(2) , xyz1(2) , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( xyz(3) , xyz1(3) , 1.e3_c_double )

  ! Check distance
  distance = geometry%distance(lonlat1(1), lonlat1(2), lonlat2(1), lonlat2(2))
  FCTEST_CHECK_CLOSE( distance , 9.971e6_c_double , 1.e3_c_double )
  distance = geometry%distance(lonlat1, lonlat2)
  FCTEST_CHECK_CLOSE( distance , 9.971e6_c_double , 1.e3_c_double )
  distance = geometry%distance(xyz1(1), xyz1(2), xyz1(3), xyz2(1), xyz2(2), xyz2(3))
  FCTEST_CHECK_CLOSE( distance , 9.971e6_c_double , 1.e3_c_double )
  distance = geometry%distance(xyz1, xyz2)
  FCTEST_CHECK_CLOSE( distance , 9.971e6_c_double , 1.e3_c_double )

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry%radius() , 3389500.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry%area() , 1.444e14_c_double , 1.e11_c_double )

  ! Finalization
  call geometry%final()

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

