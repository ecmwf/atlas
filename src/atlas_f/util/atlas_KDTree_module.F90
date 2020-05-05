! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_KDTree_module

use, intrinsic :: iso_c_binding
use fckit_object_module, only : fckit_object
use atlas_point_module
use atlas_geometry_module

implicit none

private :: fckit_object

public :: atlas_IndexKDTree

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_IndexKDTree

! Purpose :
! -------
!   *IndexKDTree* : Container of IndexKDTree

! Methods :
! -------

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_IndexKDTree__delete
  procedure :: reserve => IndexKDTree__reserve
  procedure :: insert_point => IndexKDTree__insert_point
  procedure :: insert_real => IndexKDTree__insert_real
  generic :: insert => insert_point, insert_real
  procedure :: build_only => IndexKDTree__build_only
  procedure :: build_point => IndexKDTree__build_point
  procedure :: build_real => IndexKDTree__build_real
  generic :: build => build_only, build_point, build_real
  procedure :: closestPoints_point => IndexKDTree__closestPoints_point
  procedure :: closestPoints_real => IndexKDTree__closestPoints_real
  generic :: closestPoints => closestPoints_point, closestPoints_real
  procedure :: closestPoint_point => IndexKDTree__closestPoint_point
  procedure :: closestPoint_real => IndexKDTree__closestPoint_real
  generic :: closestPoint => closestPoint_point, closestPoint_real
  procedure :: closestPointsWithinRadius_point => IndexKDTree__closestPointsWithinRadius_point
  procedure :: closestPointsWithinRadius_real => IndexKDTree__closestPointsWithinRadius_real
  generic :: closestPointsWithinRadius => closestPointsWithinRadius_point, closestPointsWithinRadius_real
  procedure :: geometry => IndexKDTree__geometry
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_IndexKDTree__final_auto
#endif

END TYPE atlas_IndexKDTree

!------------------------------------------------------------------------------

interface atlas_IndexKDTree
  module procedure atlas_IndexKDTree_cptr
  module procedure atlas_IndexKDTree__ctor
  module procedure atlas_IndexKDTree__ctor_geometry
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! IndexKDTree routines

function atlas_IndexKDTree_cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_IndexKDTree) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function atlas_IndexKDTree_cptr

function atlas_IndexKDTree__ctor() result(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  type(atlas_IndexKDTree) :: this
  call this%reset_c_ptr( atlas__IndexKDTree__new() )
end function atlas_IndexKDTree__ctor

function atlas_IndexKDTree__ctor_geometry(geometry) result(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  type(atlas_Geometry), intent(in) :: geometry
  type(atlas_IndexKDTree) :: this
  call this%reset_c_ptr( atlas__IndexKDTree__new_geometry( geometry%CPTR_PGIBUG_A ) )
end function atlas_IndexKDTree__ctor_geometry

subroutine atlas_IndexKDTree__delete(this)
  use atlas_KDTree_c_binding
  class(atlas_IndexKDTree), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__IndexKDTree__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_IndexKDTree__delete

subroutine IndexKDTree__reserve(this, size)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  integer(c_int), intent(in) :: size
  call atlas__IndexKDTree__reserve(this%CPTR_PGIBUG_A, size)
end subroutine IndexKDTree__reserve

subroutine IndexKDTree__insert_point(this, p, index)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: p
  integer(c_int), intent(in) :: index
  type(atlas_Geometry) :: geometry
  type(atlas_PointXYZ) :: p_tmp
  geometry = this%geometry()
  p_tmp = atlas_PointXYZ()
  call geometry%lonlat2xyz(p, p_tmp)
  call atlas__IndexKDTree__insert(this%CPTR_PGIBUG_A, p_tmp%CPTR_PGIBUG_A, index)
  call p_tmp%final()
end subroutine IndexKDTree__insert_point

subroutine IndexKDTree__insert_real(this, lon, lat, index)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  integer(c_int), intent(in) :: index
  type(atlas_PointLonLat) :: p
  p = atlas_PointLonLat(lon, lat)
  call this%insert(p, index)
end subroutine IndexKDTree__insert_real

subroutine IndexKDTree__build_only(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  call atlas__IndexKDTree__build(this%CPTR_PGIBUG_A)
end subroutine IndexKDTree__build_only

subroutine IndexKDTree__build_point(this, k, points, indices)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  integer(c_int), intent(in) :: k  
  type(atlas_PointLonLat), intent(in) :: points(k)
  integer(c_int), intent(in) :: indices(k)
  integer(c_int) :: i
  do i = 1, k
    call this%insert(points(i), indices(i))
  end do
  call this%build()
end subroutine IndexKDTree__build_point

subroutine IndexKDTree__build_real(this, k, lons, lats, indices)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  integer(c_int), intent(in) :: k  
  real(c_double), intent(in) :: lons(k)
  real(c_double), intent(in) :: lats(k)
  integer(c_int), intent(in) :: indices(k)
  integer(c_int) :: i
  do i = 1, k
    call this%insert(lons(i), lats(i), indices(i))
  end do
  call this%build()
end subroutine IndexKDTree__build_real

subroutine IndexKDTree__closestPoints_point(this, p, k, points, indices, distances)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: p
  integer(c_int), intent(in) :: k
  type(atlas_PointLonLat), intent(inout) :: points(k)
  integer(c_int), intent(out) :: indices(k)
  real(c_double), intent(out) :: distances(k)
  integer(c_size_t) :: k_tmp
  type(atlas_Geometry) :: geometry
  type(atlas_PointXYZ) :: p_tmp
  integer(c_int) :: i
  type(c_ptr) :: points_cptr
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  type(atlas_PointLonLat), pointer :: points_fptr(:)
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  k_tmp = int(k, c_size_t)
  geometry = this%geometry()
  p_tmp = atlas_PointXYZ()
  call geometry%lonlat2xyz(p, p_tmp)
  call atlas__IndexKDTree__closestPoints(this%CPTR_PGIBUG_A, p_tmp%CPTR_PGIBUG_A, k_tmp, points_cptr, indices_cptr, &
                                       & distances_cptr)
  call c_f_pointer(points_cptr, points_fptr, (/k/))
  call c_f_pointer(indices_cptr, indices_fptr, (/k/))
  call c_f_pointer(distances_cptr, distances_fptr, (/k/))
  indices(:) = indices_fptr(:)
  distances(:) = distances_fptr(:)
  do i = 1, k
    p_tmp = atlas_PointXYZ(c_loc(points_fptr(i)))
    points(i) = atlas_PointLonLat()
    call geometry%xyz2lonlat(p_tmp, points(i))
    call points(i)%normalise()
  end do
  call c_ptr_free(points_cptr)
  call c_ptr_free(indices_cptr)
  call c_ptr_free(distances_cptr)
  call p_tmp%final()
end subroutine IndexKDTree__closestPoints_point

subroutine IndexKDTree__closestPoints_real(this, plon, plat, k, lons, lats, indices, distances)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  integer(c_int), intent(in) :: k
  real(c_double), intent(out) :: lons(k)
  real(c_double), intent(out) :: lats(k)
  integer(c_int), intent(out) :: indices(k)
  real(c_double), intent(out) :: distances(k)
  integer(c_int) :: i
  type(atlas_PointLonLat) :: p
  type(atlas_PointLonLat) :: points(k)
  p = atlas_PointLonLat(plon, plat)
  call this%closestPoints(p, k, points, indices, distances)
  do i = 1, k
    lons(i) = points(i)%lon()
    lats(i) = points(i)%lat()
  end do
end subroutine IndexKDTree__closestPoints_real

subroutine IndexKDTree__closestPoint_point(this, p, point, index, distance)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: p
  type(atlas_PointLonLat), intent(inout) :: point
  integer(c_int), intent(out) :: index
  real(c_double), intent(out) :: distance
  type(atlas_PointLonLat) :: p_tmp
  type(c_ptr) :: point_cptr
  type(atlas_PointLonLat), pointer :: point_fptr
  p_tmp = atlas_PointLonLat(p%lon(), p%lat())
  call p_tmp%normalise()
  call atlas__IndexKDTree__closestPoint(this%CPTR_PGIBUG_A, p_tmp%CPTR_PGIBUG_A, point_cptr, index, distance)
  p_tmp = atlas_PointLonLat(point_cptr)
  point = atlas_PointLonLat(p_tmp%lon(), p_tmp%lat())
  call p_tmp%final()
end subroutine IndexKDTree__closestPoint_point

subroutine IndexKDTree__closestPoint_real(this, plon, plat, lon, lat, index, distance)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  real(c_double), intent(out) :: lon
  real(c_double), intent(out) :: lat
  integer(c_int), intent(out) :: index
  real(c_double), intent(out) :: distance
  integer(c_int) :: i
  type(atlas_PointLonLat) :: p
  type(atlas_PointLonLat) :: point
  p = atlas_PointLonLat(plon, plat)
  point = atlas_PointLonLat()
  call this%closestPoint(p, point, index, distance)
  lon = point%lon()
  lat = point%lat()
end subroutine IndexKDTree__closestPoint_real

subroutine IndexKDTree__closestPointsWithinRadius_point(this, p, radius, points, indices, distances)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: p
  real(c_double), intent(in) :: radius
  type(atlas_PointLonLat), allocatable, intent(inout) :: points(:)
  integer(c_int), allocatable, intent(inout) :: indices(:)
  real(c_double), allocatable, intent(inout) :: distances(:)
  type(atlas_PointLonLat) :: p_tmp
  integer(c_int) :: i
  integer(c_size_t) :: k
  type(c_ptr) :: points_cptr
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  type(atlas_PointLonLat), pointer :: points_fptr(:)
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  p_tmp = atlas_PointLonLat(p%lon(), p%lat())
  call p_tmp%normalise()  
  call atlas__IndexKDTree__closestPointsWithinRadius(this%CPTR_PGIBUG_A, p_tmp%CPTR_PGIBUG_A, radius, &
                                                   & points_cptr, indices_cptr, distances_cptr, k)
  call c_f_pointer(points_cptr, points_fptr, (/k/))
  call c_f_pointer(indices_cptr, indices_fptr, (/k/))
  call c_f_pointer(distances_cptr, distances_fptr, (/k/))
  if (allocated(points)) deallocate(points)
  if (allocated(indices)) deallocate(indices)
  if (allocated(distances)) deallocate(distances)
  allocate(points(k))
  allocate(indices(k))
  allocate(distances(k))
  do i = 1, int(k, c_int)
    p_tmp = atlas_PointLonLat(c_loc(points_fptr(i)))
    points(i) = atlas_PointLonLat(p_tmp%lon(), p_tmp%lat())
  end do
  indices(:) = indices_fptr(:)
  distances(:) = distances_fptr(:)
  call c_ptr_free(points_cptr)
  call c_ptr_free(indices_cptr)
  call c_ptr_free(distances_cptr)
  call p_tmp%final()
end subroutine IndexKDTree__closestPointsWithinRadius_point

subroutine IndexKDTree__closestPointsWithinRadius_real(this, plon, plat, radius, lons, lats, indices, distances)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  real(c_double), intent(in) :: radius
  real(c_double), allocatable, intent(inout) :: lons(:)
  real(c_double), allocatable, intent(inout) :: lats(:)
  integer(c_int), allocatable, intent(inout) :: indices(:)
  real(c_double), allocatable, intent(inout) :: distances(:)
  integer(c_int) :: i, k
  type(atlas_PointLonLat) :: p
  type(atlas_PointLonLat), allocatable :: points(:)
  p = atlas_PointLonLat(plon, plat)
  call this%closestPointsWithinRadius(p, radius, points, indices, distances)
  k = size(points)
  if (allocated(lons)) deallocate(lons)
  if (allocated(lats)) deallocate(lats)
  allocate(lons(k))
  allocate(lats(k))
  do i = 1, k
    lons(i) = points(i)%lon()
    lats(i) = points(i)%lat()
  end do
end subroutine IndexKDTree__closestPointsWithinRadius_real

function IndexKDTree__geometry(this) result(geometry)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  type(atlas_Geometry) :: geometry
  call geometry%reset_c_ptr( atlas__IndexKDTree__geometry( this%CPTR_PGIBUG_A ) )
end function IndexKDTree__geometry

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_IndexKDTree__final_auto(this)
  type(atlas_IndexKDTree), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_IndexKDTree__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

end module atlas_KDTree_module
