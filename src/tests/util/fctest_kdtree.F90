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
  integer(ATLAS_KIND_IDX) :: n, i, ix, iy
  integer(c_int) :: k, kk
  integer(ATLAS_KIND_IDX), allocatable :: tree_indices(:), indices(:), indices_rad(:), result_indices(:)
  real(c_double) :: lonlat(2), xyz(3), plonlat(2), pxyz(3)
  real(c_double), allocatable :: tree_lonlats(:,:), tree_distances(:)
  real(c_double), allocatable :: lonlats(:,:), distances(:)
  real(c_double), allocatable :: lons_rad(:), lats_rad(:), lonlats_rad(:,:), distances_rad(:)
  real(c_double), allocatable :: result_lonlats(:,:), result_distances(:)

  write(*,*) "test_kdtree for UnitSphere starting"

  ! Define grid
  grid = atlas_StructuredGrid("O32")

  ! Allocation
  n = grid%size()
  allocate(tree_lonlats(n,2))
  allocate(tree_indices(n))
  allocate(tree_distances(n))
  k = 6
  allocate(result_lonlats(k,2))
  allocate(result_indices(k))
  allocate(result_distances(k))
  allocate(lonlats(k,2))
  allocate(indices(k))
  allocate(distances(k))

  ! Define tree points
  i = 0
  do iy = 1, int(grid%ny(), c_int)
    do ix = 1, int(grid%nx(iy), c_int)
      i = i+1
      tree_indices(i) = i
      tree_lonlats(i,:) = grid%lonlat(ix, iy)
    end do
  end do

  ! Define test point
  plonlat = (/-76.1_c_double, -33._c_double/)

  ! Build geometry
  geometry = atlas_Geometry("UnitSphere")

  ! Define result points
  call geometry%lonlat2xyz(plonlat, pxyz)
  tree_distances = 0.0
  do i = 1, n
    call geometry%lonlat2xyz(tree_lonlats(i,:), xyz)
    tree_distances(i) = sqrt(sum((xyz-pxyz)**2))
  end do
  do i = 1, k
    result_indices(i:i) = minloc(tree_distances)
    result_distances(i) = tree_distances(result_indices(i))
    tree_distances(result_indices(i)) = huge(c_double)
    result_lonlats(i,:) = tree_lonlats(result_indices(i),:)
  end do

  ! Check interfaces with separate coordinates

  ! Check constructor with specific geometry
  kdtree = atlas_IndexKDTree(geometry)
  write(0,*) "kdtree%c_ptr() = ", c_ptr_to_loc(kdtree%CPTR_PGIBUG_A)

  ! Check reserve
  call kdtree%reserve(n)

  ! Check insert
  do i = 1, n
    call kdtree%insert(tree_lonlats(i,1), tree_lonlats(i,2), tree_indices(i))
  end do

  ! Check build only
  call kdtree%build()

  ! Check closestPoints
  call kdtree%closestPoints(plonlat(1), plonlat(2), k, indices, distances, lonlats(:,1), lonlats(:,2))
  do i = 1, k
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lonlats(i,1) , result_lonlats(i,1) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( lonlats(i,2) , result_lonlats(i,2) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-12_c_double )
  end do

  ! Check closestPoint
  call kdtree%closestPoint(plonlat(1), plonlat(2), indices(1), distances(1), lonlats(1,1), lonlats(1,2))
  FCTEST_CHECK_EQUAL( indices(1) , result_indices(1) )
  FCTEST_CHECK_CLOSE( lonlats(1,1) , result_lonlats(1,1) , 1.e-12_c_double )
  FCTEST_CHECK_CLOSE( lonlats(1,2) , result_lonlats(1,2) , 1.e-12_c_double )
  FCTEST_CHECK_CLOSE( distances(1) , result_distances(1) , 1.e-12_c_double )

  ! Check closestPoints
  call kdtree%closestPointsWithinRadius(plonlat(1), plonlat(2), 5.e-2_c_double, kk, indices_rad, distances_rad, lons_rad, lats_rad)

  FCTEST_CHECK_EQUAL( kk , 3 )
  do i = 1, kk
    FCTEST_CHECK_EQUAL( indices_rad(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lons_rad(i), result_lonlats(i,1) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( lats_rad(i) , result_lonlats(i,2) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( distances_rad(i) , result_distances(i) , 1.e-12_c_double )
  end do

  ! Finalization
  call kdtree%final()

  ! Check interfaces with vectorized coordinates

  ! Check constructor with specific geometry
  kdtree = atlas_IndexKDTree(geometry)
  write(0,*) "kdtree%c_ptr() = ", c_ptr_to_loc(kdtree%CPTR_PGIBUG_A)

  ! Check reserve
  call kdtree%reserve(n)

  ! Check insert
  do i = 1, n
    call kdtree%insert(tree_lonlats(i,:), tree_indices(i))
  end do

  ! Check build only
  call kdtree%build()

  ! Check closestPoints
  call kdtree%closestPoints(plonlat, k, indices, distances, lonlats)
  do i = 1, k
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lonlats(i,1) , result_lonlats(i,1) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( lonlats(i,2) , result_lonlats(i,2) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-12_c_double )
  end do

  ! Check closestPoint
  call kdtree%closestPoint(plonlat, indices(1), distances(1), lonlats(1,:))
  FCTEST_CHECK_EQUAL( indices(1) , result_indices(1) )
  FCTEST_CHECK_CLOSE( lonlats(1,1) , result_lonlats(1,1) , 1.e-12_c_double )
  FCTEST_CHECK_CLOSE( lonlats(1,2) , result_lonlats(1,2) , 1.e-12_c_double )
  FCTEST_CHECK_CLOSE( distances(1) , result_distances(1) , 1.e-12_c_double )

  ! Check closestPoints
  call kdtree%closestPointsWithinRadius(plonlat, 5.e-2_c_double, kk, indices_rad, distances_rad, lonlats_rad)
  FCTEST_CHECK_EQUAL( kk , 3 )
  do i = 1, kk
    FCTEST_CHECK_EQUAL( indices_rad(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lonlats_rad(i,1) , result_lonlats(i,1) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( lonlats_rad(i,2) , result_lonlats(i,2) , 1.e-12_c_double )
    FCTEST_CHECK_CLOSE( distances_rad(i) , result_distances(i) , 1.e-12_c_double )
  end do

  ! Finalization
  call kdtree%final()
  call geometry%final()

  ! Check constructor without geometry (Earth)
  kdtree = atlas_IndexKDTree()
  write(0,*) "kdtree%c_ptr() = ", c_ptr_to_loc(kdtree%CPTR_PGIBUG_A)

  ! Check geometry accessor
  geometry = kdtree%geometry()
  FCTEST_CHECK_EQUAL( geometry%radius() , 6371229._c_double )

  ! Define result points
  call geometry%lonlat2xyz(plonlat(1), plonlat(2), pxyz(1), pxyz(2), pxyz(3))
  tree_distances = 0.0
  do i = 1, n
    call geometry%lonlat2xyz(tree_lonlats(i,1), tree_lonlats(i,2), xyz(1), xyz(2), xyz(3))
    tree_distances(i) = sqrt(sum((xyz-pxyz)**2))
  end do
  do i = 1, k
    result_indices(i:i) = minloc(tree_distances)
    result_distances(i) = tree_distances(result_indices(i))
    tree_distances(result_indices(i)) = huge(c_double)
    result_lonlats(i,:) = tree_lonlats(result_indices(i),:)
  end do

  ! Check interfaces with separate coordinates

  ! Check build with list
  call kdtree%build(n, tree_lonlats(:,1), tree_lonlats(:,2), tree_indices)

  ! Check closestPoints
  call kdtree%closestPoints(plonlat(1), plonlat(2), k, indices, distances, lonlats(:,1), lonlats(:,2))
  do i = 1, k
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lonlats(i,1) , result_lonlats(i,1) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( lonlats(i,2) , result_lonlats(i,2) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-10_c_double )
  end do

  ! Check closestPoint
  call kdtree%closestPoint(plonlat(1), plonlat(2), indices(1), distances(1), lonlats(1,1), lonlats(1,2))
  FCTEST_CHECK_EQUAL( indices(1) , result_indices(1) )
  FCTEST_CHECK_CLOSE( lonlats(1,1) , result_lonlats(1,1) , 1.e-10_c_double )
  FCTEST_CHECK_CLOSE( lonlats(1,2) , result_lonlats(1,2) , 1.e-10_c_double )
  FCTEST_CHECK_CLOSE( distances(1) , result_distances(1) , 1.e-10_c_double )

  ! Check closestPoints
  call kdtree%closestPointsWithinRadius(plonlat(1), plonlat(2), 3.5e5_c_double, kk, indices_rad, distances_rad, lons_rad, lats_rad)
  FCTEST_CHECK_EQUAL( kk , 4 )
  do i = 1, kk
    FCTEST_CHECK_EQUAL( indices_rad(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lons_rad(i) , result_lonlats(i,1) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( lats_rad(i) , result_lonlats(i,2) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( distances_rad(i) , result_distances(i) , 1.e-10_c_double )
  end do

  ! Finalization
  call kdtree%final()
  call geometry%final()

  ! Check interfaces with vectorized coordinates

  ! Check constructor without geometry (Earth)
  kdtree = atlas_IndexKDTree()
  write(0,*) "kdtree%c_ptr() = ", c_ptr_to_loc(kdtree%CPTR_PGIBUG_A)

  ! Check build with list
  call kdtree%build(n, tree_lonlats, tree_indices)

  ! Check closestPoints
  call kdtree%closestPoints(plonlat, k, indices, distances, lonlats)
  do i = 1, k
    FCTEST_CHECK_EQUAL( indices(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lonlats(i,1) , result_lonlats(i,1) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( lonlats(i,2) , result_lonlats(i,2) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( distances(i) , result_distances(i) , 1.e-10_c_double )
  end do

  ! Check closestPoint
  call kdtree%closestPoint(plonlat, indices(1), distances(1), lonlats(1,:))
  FCTEST_CHECK_EQUAL( indices(1) , result_indices(1) )
  FCTEST_CHECK_CLOSE( lonlats(1,1) , result_lonlats(1,1) , 1.e-10_c_double )
  FCTEST_CHECK_CLOSE( lonlats(1,2) , result_lonlats(1,2) , 1.e-10_c_double )
  FCTEST_CHECK_CLOSE( distances(1) , result_distances(1) , 1.e-10_c_double )

  ! Check closestPoints
  call kdtree%closestPointsWithinRadius(plonlat, 3.5e5_c_double, kk, indices_rad, distances_rad, lonlats_rad)
  FCTEST_CHECK_EQUAL( kk , 4 )
  do i = 1, kk
    FCTEST_CHECK_EQUAL( indices_rad(i) , result_indices(i) )
    FCTEST_CHECK_CLOSE( lonlats_rad(i,1) , result_lonlats(i,1) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( lonlats_rad(i,2) , result_lonlats(i,2) , 1.e-10_c_double )
    FCTEST_CHECK_CLOSE( distances_rad(i) , result_distances(i) , 1.e-10_c_double )
  end do

  ! Finalization
  call kdtree%final()

  ! Release memory
  call grid%final()
  deallocate(tree_lonlats)
  deallocate(tree_indices)
  deallocate(tree_distances)
  deallocate(lonlats)
  deallocate(indices)
  deallocate(distances)
  deallocate(lons_rad)
  deallocate(lats_rad)
  deallocate(lonlats_rad)
  deallocate(indices_rad)
  deallocate(distances_rad)
  deallocate(result_lonlats)
  deallocate(result_indices)
  deallocate(result_distances)

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

