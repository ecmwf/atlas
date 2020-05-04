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

  type(atlas_Geometry) :: geometry, geometry_ptr
  type(atlas_PointLonLat) :: p1LonLat, p2LonLat, p1LonLat_test
  type(atlas_PointXYZ) :: p1XYZ, p2XYZ, p1XYZ_test
  real(c_double) :: lon, lat, x, y, z, d

  write(*,*) "test_geometry for UnitSphere starting"

  ! Check constructor
  geometry = atlas_Geometry("UnitSphere")
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Define points
  p1LonLat = atlas_PointLonLat( -71.6_c_double, -33._c_double)
  p2LonLat = atlas_PointLonLat( 121.8_c_double, 31.4_c_double)
  p1XYZ = atlas_PointXYZ( 2.647e-1_c_double, -7.958e-1_c_double, -5.446e-1_c_double)
  p2XYZ = atlas_PointXYZ( -4.498e-1_c_double, 7.254e-1_c_double, 5.210e-1_c_double)

  ! Check conversion from spherical to cartesian
  p1XYZ_test = atlas_PointXYZ()
  call geometry%lonlat2xyz(p1LonLat, p1XYZ_test)
  FCTEST_CHECK_CLOSE( p1XYZ_test%x() , p1XYZ%x(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( p1XYZ_test%y() , p1XYZ%y(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( p1XYZ_test%z() , p1XYZ%z(), 1.e3_c_double )

  ! Check conversion from cartesian to spherical
  p1LonLat_test = atlas_PointLonLat()
  call geometry%xyz2lonlat(p1XYZ, p1LonLat_test)
  FCTEST_CHECK_CLOSE( p1LonLat_test%lon() , p1LonLat%lon(), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( p1LonLat_test%lat() , p1LonLat%lat(), 1.e-2_c_double )

  ! Check distance
  FCTEST_CHECK_CLOSE( geometry%distance(p1LonLat, p2LonLat) , 2.942_c_double , 1.e-3_c_double )
  FCTEST_CHECK_CLOSE( geometry%distance(p1XYZ, p2XYZ) , 2.942_c_double , 1.e-3_c_double )

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry%radius() , 1.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry%area() , 1.257e1_c_double , 1.e-2_c_double )

  ! Finalization
  call geometry%final()
  call p1LonLat%final()
  call p2LonLat%final()
  call p1XYZ%final()
  call p2XYZ%final()

  write(*,*) "test_geometry for Earth starting"

  ! Check constructor for Earth
  geometry = atlas_Geometry("Earth")
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Define points
  p1LonLat = atlas_PointLonLat( -71.6_c_double, -33._c_double)
  p2LonLat = atlas_PointLonLat( 121.8_c_double, 31.4_c_double)
  p1XYZ = atlas_PointXYZ( 1.687e6_c_double, -5.070e6_c_double, -3.470e6_c_double)
  p2XYZ = atlas_PointXYZ( -2.866e6_c_double, 4.622e6_c_double, 3.319e6_c_double)

  ! Check conversion from cartesian to spherical
  p1LonLat_test = atlas_PointLonLat()
  call geometry%xyz2lonlat(p1XYZ, p1LonLat_test)
  FCTEST_CHECK_CLOSE( p1LonLat_test%lon() , p1LonLat%lon(), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( p1LonLat_test%lat() , p1LonLat%lat(), 1.e-2_c_double )

  ! Check conversion from spherical to cartesian
  p1XYZ_test = atlas_PointXYZ()
  call geometry%lonlat2xyz(p1LonLat, p1XYZ_test)
  FCTEST_CHECK_CLOSE( p1XYZ_test%x() , p1XYZ%x(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( p1XYZ_test%y() , p1XYZ%y(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( p1XYZ_test%z() , p1XYZ%z(), 1.e3_c_double )

  ! Check distance
  FCTEST_CHECK_CLOSE( geometry%distance(p1LonLat, p2LonLat) , 1.874e7_c_double , 1.e4_c_double )
  FCTEST_CHECK_CLOSE( geometry%distance(p1XYZ, p2XYZ) , 1.874e7_c_double , 1.e4_c_double )

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry%radius() , 6371229.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry%area() , 5.101e14_c_double , 1.e11_c_double )

  ! Finalization
  call geometry%final()
  call p1LonLat%final()
  call p2LonLat%final()
  call p1XYZ%final()
  call p2XYZ%final()

  write(*,*) "test_geometry for another planet (Mars here) starting"

  ! Check constructor for another planet (Mars here)
  geometry = atlas_Geometry( 3389500.0_c_double )
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Define points
  p1LonLat = atlas_PointLonLat( -71.6_c_double, -33._c_double)
  p2LonLat = atlas_PointLonLat( 121.8_c_double, 31.4_c_double)  
  p1XYZ = atlas_PointXYZ( 8.973e5_c_double, -2.697e6_c_double, -1.846e6_c_double)
  p2XYZ = atlas_PointXYZ( -1.525e6_c_double, 2.459e6_c_double, 1.766e6_c_double)

  ! Check conversion from cartesian to spherical
  p1LonLat_test = atlas_PointLonLat()
  call geometry%xyz2lonlat(p1XYZ, p1LonLat_test)
  FCTEST_CHECK_CLOSE( p1LonLat_test%lon() , p1LonLat%lon(), 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( p1LonLat_test%lat() , p1LonLat%lat(), 1.e-2_c_double )

  ! Check conversion from spherical to cartesian
  p1XYZ_test = atlas_PointXYZ()
  call geometry%lonlat2xyz(p1LonLat, p1XYZ_test)
  FCTEST_CHECK_CLOSE( p1XYZ_test%x() , p1XYZ%x(), 1.e2_c_double )
  FCTEST_CHECK_CLOSE( p1XYZ_test%y() , p1XYZ%y(), 1.e3_c_double )
  FCTEST_CHECK_CLOSE( p1XYZ_test%z() , p1XYZ%z(), 1.e3_c_double )

  ! Check distance
  FCTEST_CHECK_CLOSE( geometry%distance(p1LonLat, p2LonLat) , 9.971e6_c_double , 1.e3_c_double )
  FCTEST_CHECK_CLOSE( geometry%distance(p1XYZ, p2XYZ) , 9.971e6_c_double , 1.e3_c_double )

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry%radius() , 3389500.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry%area() , 1.444e14_c_double , 1.e11_c_double )

  ! Finalization
  call geometry%final()
  call p1LonLat%final()
  call p2LonLat%final()
  call p1XYZ%final()
  call p2XYZ%final()

  write(*,*) "test_geometry for real interfaces"

  ! Check constructor
  geometry = atlas_Geometry("UnitSphere")
  write(0,*) "geometry%c_ptr() = ", c_ptr_to_loc(geometry%CPTR_PGIBUG_A)

  ! Check conversion from spherical to cartesian
  call geometry%lonlat2xyz(-71.6_c_double, -33._c_double, x, y, z)
  FCTEST_CHECK_CLOSE( x , 2.647e-1_c_double, 1.e3_c_double )
  FCTEST_CHECK_CLOSE( y , -7.958e-1_c_double, 1.e3_c_double )
  FCTEST_CHECK_CLOSE( z , -5.446e-1_c_double, 1.e3_c_double )

  ! Check conversion from cartesian to spherical
  call geometry%xyz2lonlat(2.647e-1_c_double, -7.958e-1_c_double, -5.446e-1_c_double, lon, lat)
  FCTEST_CHECK_CLOSE( lon , -71.6_c_double, 1.e-2_c_double )
  FCTEST_CHECK_CLOSE( lat , -33._c_double, 1.e-2_c_double )

  ! Check distance
  d = geometry%distance(-71.6_c_double, -33._c_double, 121.8_c_double, 31.4_c_double)
  FCTEST_CHECK_CLOSE( d , 2.942_c_double , 1.e-3_c_double )
  d = geometry%distance(2.647e-1_c_double, -7.958e-1_c_double, -5.446e-1_c_double, &
                      & -4.498e-1_c_double, 7.254e-1_c_double, 5.210e-1_c_double)
  FCTEST_CHECK_CLOSE( d , 2.942_c_double , 1.e-3_c_double )

  ! Check constructor with pointer
  geometry_ptr = atlas_Geometry(geometry%c_ptr())
  write(0,*) "geometry_ptr%c_ptr() = ", c_ptr_to_loc(geometry_ptr%CPTR_PGIBUG_A)

  ! Check radius
  FCTEST_CHECK_EQUAL( geometry_ptr%radius() , 1.0_c_double )

  ! Check area
  FCTEST_CHECK_CLOSE( geometry_ptr%area() , 1.257e1_c_double , 1.e-2_c_double )

  ! Finalization
  call geometry_ptr%final()
  call geometry%final()
  call p1LonLat%final()
  call p2LonLat%final()
  call p1XYZ%final()
  call p2XYZ%final()

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

