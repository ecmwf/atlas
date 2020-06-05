! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! @author Willem Deconinck

#include "fckit/fctest.h"

module fixture_reg_grid
use atlas_module
use, intrinsic :: iso_c_binding, only : c_double
implicit none
contains

subroutine print_spec( grid )
   type(atlas_Grid) :: grid
   type(atlas_Config) :: spec
   spec = grid%spec()
   write(0,*) spec%json()
   call spec%final()
end subroutine

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_regional_Grid, fixture_reg_grid)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

TEST( simple_regional_grid )
#if 1
    type(atlas_Grid) :: grid1, grid2
    grid1 = atlas_RegionalGrid( nx=11, ny=11, xy_min=[0._dp,-5._dp], xy_max=[10._dp,5._dp] )
    grid2 = atlas_RegionalGrid( nx=11, ny=11, north=5._dp, west=0._dp, south=-5._dp, east=10._dp )
    call print_spec(grid1)
    call print_spec(grid2)
    FCTEST_CHECK_EQUAL( grid1%uid() , "0aae0dfb8141e8c70ef63108e6c59ac1" )
    FCTEST_CHECK_EQUAL( grid1%uid() , grid2%uid() )
    call grid1%final()
    call grid2%final()
#endif
END_TEST

TEST( rotated_regional_grid )
#if 1
    type(atlas_Grid) :: grid1, grid2
    grid1 = atlas_RegionalGrid(  nx=11, ny=11, xy_min=[0._dp,-5._dp], xy_max=[10._dp,5._dp], &
               & projection=atlas_RotatedLonLatProjection([40._dp,20._dp]) )
    grid2 = atlas_RegionalGrid( nx=11, ny=11, north=5._dp, west=0._dp, south=-5._dp, east=10._dp, &
               & projection=atlas_RotatedLonLatProjection([40._dp,20._dp]) )
    call print_spec(grid1)
    FCTEST_CHECK_EQUAL( grid1%uid() , grid2%uid() )
    call grid1%final()
    call grid2%final()
#endif
END_TEST

TEST( lambert_grid )
#if 1
    type(atlas_Grid) :: grid
    grid = atlas_RegionalGrid( nx=11, ny=11, dx=10000._dp,dy=10000._dp, &
               & xy_min=[-50000._dp,-50000._dp],  &
               & projection=atlas_LambertConformalConicProjection(4.0_dp,50._dp) )
    call print_spec(grid)
    call grid%final()
#endif
END_TEST

TEST( test_regional_lonlat_grid_MF )
#if 1
  ! Grid provided by Philippe Marguinaud, Meteo France
  type(atlas_StructuredGrid) :: grid
  real(c_double), parameter :: zlonw = -20., zlone = +10., zlats = 10., zlatn = 50.
  integer(c_int), parameter :: ilons = 30, ilats = 40
  real(c_double), parameter :: tol = 1.e-5_dp

  grid = atlas_RegionalGrid( nx=ilons, ny=ilats, north=zlatn, west=zlonw, south=zlats, east=zlone )

  FCTEST_CHECK_EQUAL( grid%size(), ilons*ilats )

  FCTEST_CHECK_CLOSE( grid%lonlat(1,1),         ([-0.2000000000E+02_dp, 0.1000000000E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(2,1),         ([-0.1896551724E+02_dp, 0.1000000000E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(ilons,1),     ([ 0.1000000000E+02_dp, 0.1000000000E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(1,2),         ([-0.2000000000E+02_dp, 0.1102564103E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(ilons,ilats), ([ 0.1000000000E+02_dp, 0.5000000000E+02_dp]), tol)
#endif
END_TEST



TEST( test_regional_lambert_grid_MF )
#if 1
  ! Grid provided by Philippe Marguinaud, Meteo France
  type(atlas_StructuredGrid) :: grid

  integer(c_int), parameter :: ndlon=64
  integer(c_int), parameter :: ndglg=64
  integer(c_int), parameter :: nux=53
  integer(c_int), parameter :: nuy=53
  real(c_double), parameter :: dxinmetres=50000.
  real(c_double), parameter :: dyinmetres=50000.
  real(c_double), parameter :: xmin = int(real(-nux,c_double)/2.) * dxinmetres
  real(c_double), parameter :: ymin = int(real(-nuy,c_double)/2.) * dyinmetres
  real(c_double), parameter :: ladindegrees=46.2
  real(c_double), parameter :: latin1indegrees=46.2
  real(c_double), parameter :: latin2indegrees=46.2
  real(c_double), parameter :: lovindegrees=2.0
  real(c_double), parameter :: tol = 1.e-5_dp

  grid = atlas_RegionalGrid( nx=ndlon, ny=ndglg, xy_min=[xmin,ymin], &
               & dx=dxinmetres, dy=dyinmetres, &
               & projection=atlas_LambertConformalConicProjection(lovindegrees,ladindegrees, &
               &            latin1indegrees,latin2indegrees) )

  FCTEST_CHECK_EQUAL( grid%size(), ndglg*ndlon )

  FCTEST_CHECK_CLOSE( grid%lonlat(1,1),         ([-0.1178699484E+02_dp, 0.3358923512E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(2,1),         ([-0.1126673471E+02_dp, 0.3366377557E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(ndlon,1),     ([ 0.2142256121E+02_dp, 0.3258651702E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(1,2),         ([-0.1187876836E+02_dp, 0.3402241700E+02_dp]), tol)
  FCTEST_CHECK_CLOSE( grid%lonlat(ndlon,ndglg), ([ 0.3452466369E+02_dp, 0.5925747619E+02_dp]), tol)

  call grid%final()
#endif
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

