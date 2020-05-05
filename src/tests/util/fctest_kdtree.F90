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

  type(atlas_Geometry) :: geometry, geometry_test
  type(atlas_IndexKDTree) :: kdtree_point, kdtree_real
  type(atlas_StructuredGrid) :: grid
  integer(c_int) :: n, i, ix, iy, k
  integer(c_int), allocatable :: tree_indices(:), indices(:), result_indices(:)
  real(c_double) :: lonlat(2), plon, plat
  real(c_double), allocatable :: tree_lons(:), tree_lats(:), tree_distances(:)
  real(c_double), allocatable :: lons(:), lats(:), distances(:)
  real(c_double), allocatable :: result_lons(:), result_lats(:), result_distances(:)
  type(atlas_PointLonLat) :: p
  type(atlas_PointLonLat), allocatable :: tree_points(:), points(:), result_points(:)

  write(*,*) "test_kdtree for UnitSphere starting"

  ! Define grid and geometry
  grid = atlas_StructuredGrid("O32")
  geometry = atlas_Geometry("UnitSphere")

  ! Allocation
  n = grid%size()
  allocate(tree_indices(n))
  allocate(tree_lons(n))
  allocate(tree_lats(n))
  allocate(tree_points(n))
  allocate(tree_distances(n))
  k = 6
  allocate(result_points(k))
  allocate(result_lons(k))
  allocate(result_lats(k))
  allocate(result_indices(k))
  allocate(result_distances(k))
  allocate(points(k))
  allocate(indices(k))
  allocate(distances(k))
  allocate(lons(k))
  allocate(lats(k))

  ! Define tree points
  i = 0
  do iy = 1, int(grid%ny(), c_int)
    do ix = 1, int(grid%nx(iy), c_int)
      i = i+1
      tree_indices(i) = i
      lonlat = grid%lonlat(ix, iy)
      tree_lons(i) = lonlat(1)
      tree_lats(i) = lonlat(2)
      tree_points(i) = atlas_PointLonLat(tree_lons(i), tree_lats(i))
    end do
  end do

  ! Define test point
  plon = -76.1_c_double
  plat = -33._c_double
  p = atlas_PointLonLat(plon, plat)

  ! Define result points
  tree_distances = 0.0
  do i = 1, n
    tree_distances(i) = geometry%distance(p, tree_points(i))
  end do
  do i = 1, k
    result_indices(i:i) = minloc(tree_distances)
    result_distances(i) = tree_distances(result_indices(i))
    tree_distances(result_indices(i)) = huge(c_double)
    result_lons(i) = tree_lons(result_indices(i))
    result_lats(i) = tree_lats(result_indices(i))
    result_points(i) = atlas_PointLonLat(result_lons(i), result_lats(i))
  end do

  ! Check constructor
  kdtree_point = atlas_IndexKDTree(geometry)
  kdtree_real = atlas_IndexKDTree(geometry)
  write(0,*) "kdtree_point%c_ptr() = ", c_ptr_to_loc(kdtree_point%CPTR_PGIBUG_A)
  write(0,*) "kdtree_real%c_ptr() = ", c_ptr_to_loc(kdtree_real%CPTR_PGIBUG_A)

  ! Check geometry accessor
  geometry_test = kdtree_point%geometry()
  FCTEST_CHECK_EQUAL( geometry_test%radius() , 1.0_c_double )

  ! Check reserve
  call kdtree_point%reserve(n)

  ! Check insert
  do i = 1, n
    call kdtree_point%insert(tree_points(i), tree_indices(i))
    call kdtree_real%insert(tree_lons(i), tree_lats(i), tree_indices(i))
  end do

  ! Check build only
  call kdtree_point%build()
  call kdtree_real%build()

  ! Check closestPoints
  call kdtree_point%closestPoints(p, k, points, indices, distances)
  do i = 1, k
    FCTEST_CHECK_CLOSE( points(i)%lon() , result_lons(i), 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( points(i)%lat() , result_lats(i), 1.e-12_c_double )
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-3_c_double )
  end do
  call kdtree_real%closestPoints(plon, plat, k, lons, lats, indices, distances)
  do i = 1, k
    FCTEST_CHECK_CLOSE( points(i)%lon() , result_lons(i), 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( points(i)%lat() , result_lats(i), 1.e-12_c_double )
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-3_c_double )
  end do
stop
  ! Check closestPoint
  call kdtree_point%closestPoint(p, points(1), indices(1), distances(1))
  FCTEST_CHECK_EQUAL( points(1)%lon() , result_lons(1) )
  FCTEST_CHECK_EQUAL( points(1)%lat() , result_lats(1) )
  FCTEST_CHECK_EQUAL( indices(1) , result_indices(1) )
  FCTEST_CHECK_CLOSE( distances(1) , result_distances(1) , 1.e-3_c_double )
  call kdtree_real%closestPoint(plon, plat, lons(1), lats(1), indices(1), distances(1))
  FCTEST_CHECK_EQUAL( lons(1) , result_lons(1) )
  FCTEST_CHECK_EQUAL( lats(1) , result_lats(1) )
  FCTEST_CHECK_EQUAL( indices(1) , result_indices(1) )
  FCTEST_CHECK_CLOSE( distances(1) , result_distances(1) , 1.e-3_c_double )

  ! Check closestPoints
  call kdtree_point%closestPointsWithinRadius(p, 3.5_c_double, points, indices, distances)
  k = size(points)
  FCTEST_CHECK_EQUAL( k , 3 )
  do i = 1, k
    FCTEST_CHECK_EQUAL( points(i)%lon() , result_lons(i) )
    FCTEST_CHECK_EQUAL( points(i)%lat() , result_lats(i) )
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-3_c_double )
  end do
  call kdtree_real%closestPointsWithinRadius(plon, plat, 3.5_c_double, lons, lats, indices, distances)
  k = size(points)
  FCTEST_CHECK_EQUAL( k , 3 )
  do i = 1, k
    FCTEST_CHECK_EQUAL( lons(i) , result_lons(i) )
    FCTEST_CHECK_EQUAL( lats(i) , result_lats(i) )
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-3_c_double )
  end do

  ! Finalization
  call geometry%final()
  call kdtree_point%final()
  call kdtree_real%final()

  ! Check constructor without geometry (Earth)
  kdtree_point = atlas_IndexKDTree()
  kdtree_real = atlas_IndexKDTree()
  write(0,*) "kdtree_point%c_ptr() = ", c_ptr_to_loc(kdtree_point%CPTR_PGIBUG_A)
  write(0,*) "kdtree_real%c_ptr() = ", c_ptr_to_loc(kdtree_real%CPTR_PGIBUG_A)

  ! Check geometry accessor
  geometry_test = kdtree_point%geometry()
  FCTEST_CHECK_EQUAL( geometry_test%radius() , 6371229._c_double )

  ! Check build with arrays
  call kdtree_point%build(n, tree_points, tree_indices)
  call kdtree_real%build(n, tree_lons, tree_lats, tree_indices)

  ! Check closestPoints
  call kdtree_point%closestPoints(p, k, points, indices, distances)
  do i = 1, k
    FCTEST_CHECK_EQUAL( points(i)%lon() , result_lons(i) )
    FCTEST_CHECK_EQUAL( points(i)%lat() , result_lats(i) )
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-3_c_double )
  end do
  call kdtree_real%closestPoints(plon, plat, k, lons, lats, indices, distances)
  do i = 1, k
    FCTEST_CHECK_EQUAL( lons(i) , result_lons(i) )
    FCTEST_CHECK_EQUAL( lats(i) , result_lats(i) )
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-3_c_double )
  end do

  ! Finalization
  call geometry%final()
  call geometry_test%final()
  call grid%final()
  call kdtree_point%final()
  call kdtree_real%final()
  deallocate(tree_indices)
  deallocate(indices)
  deallocate(result_indices)
  deallocate(tree_lons)
  deallocate(tree_lats)
  deallocate(lons)
  deallocate(lats)
  deallocate(distances)
  deallocate(result_lons)
  deallocate(result_lats)
  deallocate(result_distances)
  call p%final()
  do i = 1, size(tree_points)
    call tree_points(i)%final()
  end do
  deallocate(tree_points)
  do i = 1, size(points)
    call points(i)%final()
  end do
  deallocate(points)
  do i = 1, size(result_points)
    call result_points(i)%final()
  end do
  deallocate(result_points)

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

