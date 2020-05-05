! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the KDTree Datastructure
! @author Benjamin Menetrier

#include "fckit/fctest.h"
#include "atlas/atlas_f.h"

! -----------------------------------------------------------------------------

module fcta_KDTree_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_KDTree_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_KDTree,fcta_KDTree_fixture)

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

TEST( test_kdtree )
use fckit_log_module
use fckit_c_interop_module
implicit none

  type(atlas_Geometry) :: geometry
  type(atlas_IndexKDTree) :: kdtree
  type(atlas_StructuredGrid) :: grid
  integer(c_int) :: n, i, ix, iy, k, kk
  integer(c_int), allocatable :: tree_indices(:), indices(:), indices_rad(:), result_indices(:)
  real(c_double) :: lonlat(2), xyz(3), plon, plat, pxyz(3)
  real(c_double), allocatable :: tree_lons(:), tree_lats(:), tree_distances(:)
  real(c_double), allocatable :: lons(:), lats(:), distances(:)
  real(c_double), allocatable :: lons_rad(:), lats_rad(:), distances_rad(:)
  real(c_double), allocatable :: result_lons(:), result_lats(:), result_distances(:)

  write(*,*) "test_kdtree for UnitSphere starting"

  ! Define grid
  grid = atlas_StructuredGrid("O32")

  ! Allocation
  n = grid%size()
  allocate(tree_lons(n))
  allocate(tree_lats(n))
  allocate(tree_indices(n))
  allocate(tree_distances(n))
  k = 6
  allocate(result_lons(k))
  allocate(result_lats(k))
  allocate(result_indices(k))
  allocate(result_distances(k))
  allocate(lons(k))
  allocate(lats(k))
  allocate(indices(k))
  allocate(distances(k))

  ! Define tree points
  i = 0
  do iy = 1, int(grid%ny(), c_int)
    do ix = 1, int(grid%nx(iy), c_int)
      i = i+1
      tree_indices(i) = i
      lonlat = grid%lonlat(ix, iy)
      tree_lons(i) = lonlat(1)
      tree_lats(i) = lonlat(2)
    end do
  end do

  ! Define test point
  plon = -76.1_c_double
  plat = -33._c_double

  ! Check constructor with specific geometry
  geometry = atlas_Geometry("UnitSphere")
  kdtree = atlas_IndexKDTree(geometry)
!  write(0,*) "kdtree%c_ptr() = ", c_ptr_to_loc(kdtree%CPTR_PGIBUG_A)

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

