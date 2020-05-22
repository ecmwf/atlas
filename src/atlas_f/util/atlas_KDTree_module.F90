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
use fckit_owned_object_module, only : fckit_owned_object
use atlas_geometry_module

implicit none

private :: fckit_owned_object

public :: atlas_IndexKDTree

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_IndexKDTree

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
  procedure :: insert_separate_coords => IndexKDTree__insert_separate_coords
  procedure :: insert_vectorized_coords => IndexKDTree__insert_vectorized_coords
  generic :: insert => insert_separate_coords, insert_vectorized_coords
  procedure :: build_only => IndexKDTree__build_only
  procedure :: build_list_separate_coords => IndexKDTree__build_list_separate_coords
  procedure :: build_list_vectorized_coords => IndexKDTree__build_list_vectorized_coords
  generic :: build => build_only, build_list_separate_coords, build_list_vectorized_coords
  procedure :: closestPoints_separate_coords => IndexKDTree__closestPoints_separate_coords
  procedure :: closestPoints_vectorized_coords => IndexKDTree__closestPoints_vectorized_coords
  generic :: closestPoints => closestPoints_separate_coords, closestPoints_vectorized_coords
  procedure :: closestPoint_separate_coords => IndexKDTree__closestPoint_separate_coords
  procedure :: closestPoint_vectorized_coords => IndexKDTree__closestPoint_vectorized_coords
  generic :: closestPoint => closestPoint_separate_coords, closestPoint_vectorized_coords
  procedure :: closestPointsWithinRadius_separate_coords => IndexKDTree__closestPointsWithinRadius_separate_coords
  procedure :: closestPointsWithinRadius_vectorized_coords => IndexKDTree__closestPointsWithinRadius_vectorized_coords
  generic :: closestPointsWithinRadius => closestPointsWithinRadius_separate_coords, closestPointsWithinRadius_vectorized_coords
  procedure :: geometry => IndexKDTree__geometry
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_IndexKDTree__final_auto
#endif

END TYPE atlas_IndexKDTree

!------------------------------------------------------------------------------

interface atlas_IndexKDTree
  module procedure atlas_IndexKDTree__ctor
  module procedure atlas_IndexKDTree__ctor_geometry
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! IndexKDTree routines

function atlas_IndexKDTree__ctor() result(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  type(atlas_IndexKDTree) :: this
  call this%reset_c_ptr( atlas__IndexKDTree__new() )
  call this%return()
end function atlas_IndexKDTree__ctor

function atlas_IndexKDTree__ctor_geometry(geometry) result(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  type(atlas_Geometry), intent(in) :: geometry
  type(atlas_IndexKDTree) :: this
  call this%reset_c_ptr( atlas__IndexKDTree__new_geometry( geometry%CPTR_PGIBUG_A ) )
  call this%return()
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

subroutine IndexKDTree__insert_separate_coords(this, lon, lat, index)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  integer(c_int), intent(in) :: index
  call atlas__IndexKDTree__insert(this%CPTR_PGIBUG_A, lon, lat, index)
end subroutine IndexKDTree__insert_separate_coords

subroutine IndexKDTree__insert_vectorized_coords(this, lonlat, index)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: lonlat(2)
  integer(c_int), intent(in) :: index
  call atlas__IndexKDTree__insert(this%CPTR_PGIBUG_A, lonlat(1), lonlat(2), index)
end subroutine IndexKDTree__insert_vectorized_coords

subroutine IndexKDTree__build_only(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  call atlas__IndexKDTree__build(this%CPTR_PGIBUG_A)
end subroutine IndexKDTree__build_only

subroutine IndexKDTree__build_list_separate_coords(this, k, lons, lats, indices)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  integer(c_int), intent(in) :: k
  real(c_double), intent(in) :: lons(k)
  real(c_double), intent(in) :: lats(k)
  integer(c_int), intent(in), optional :: indices(k)
  integer(c_int) :: i, index
  do i = 1, k
    if (present(indices)) then
      index = indices(i)
    else
      index = i
    end if
    call this%insert(lons(i), lats(i), index)
  end do
  call this%build()
end subroutine IndexKDTree__build_list_separate_coords

subroutine IndexKDTree__build_list_vectorized_coords(this, k, lonlats, indices)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  integer(c_int), intent(in) :: k  
  real(c_double), intent(in) :: lonlats(k,2)
  integer(c_int), intent(in), optional :: indices(k)
  integer(c_int) :: i, index
  do i = 1, k
    if (present(indices)) then
      index = indices(i)
    else
      index = i
    end if
    call this%insert(lonlats(i,:), index)
  end do
  call this%build()
end subroutine IndexKDTree__build_list_vectorized_coords

subroutine IndexKDTree__closestPoints_separate_coords(this, plon, plat, k, indices, distances, lons, lats)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  integer(c_int), intent(in) :: k
  integer(c_int), intent(out) :: indices(k)
  real(c_double), intent(out), optional :: distances(k)
  real(c_double), intent(out), optional :: lons(k)
  real(c_double), intent(out), optional :: lats(k)
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  type(c_ptr) :: lons_cptr
  type(c_ptr) :: lats_cptr
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  real(c_double), pointer :: lons_fptr(:)
  real(c_double), pointer :: lats_fptr(:)
  if ( k > 0 ) then
    call atlas__IndexKDTree__closestPoints(this%CPTR_PGIBUG_A, plon, plat, int(k, c_size_t), lons_cptr, lats_cptr, indices_cptr, &
                                         & distances_cptr)
    call c_f_pointer(indices_cptr, indices_fptr, (/k/))
    indices(:) = indices_fptr(:)
    if (present(distances)) then
      call c_f_pointer(distances_cptr, distances_fptr, (/k/))
      distances(:) = distances_fptr(:)
    end if
    if (present(lons)) then
      call c_f_pointer(lons_cptr, lons_fptr, (/k/))
      lons(:) = lons_fptr(:)
    end if
    if (present(lats)) then
      call c_f_pointer(lats_cptr, lats_fptr, (/k/))
      lats(:) = lats_fptr(:)
    end if
    call c_ptr_free(lons_cptr)
    call c_ptr_free(lats_cptr)
    call c_ptr_free(indices_cptr)
    call c_ptr_free(distances_cptr)
  end if
end subroutine IndexKDTree__closestPoints_separate_coords

subroutine IndexKDTree__closestPoints_vectorized_coords(this, point, k, indices, distances, lonlats)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: point(2)
  integer(c_int), intent(in) :: k
  integer(c_int), intent(out) :: indices(k)
  real(c_double), intent(out), optional :: distances(k)
  real(c_double), intent(out), optional :: lonlats(k,2)
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  type(c_ptr) :: lons_cptr
  type(c_ptr) :: lats_cptr
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  real(c_double), pointer :: lons_fptr(:)
  real(c_double), pointer :: lats_fptr(:)
  if ( k > 0 ) then
    call atlas__IndexKDTree__closestPoints(this%CPTR_PGIBUG_A, point(1), point(2), int(k, c_size_t), lons_cptr, lats_cptr, &
                                         & indices_cptr, distances_cptr)
    call c_f_pointer(indices_cptr, indices_fptr, (/k/))
    indices(:) = indices_fptr(:)
    if (present(distances)) then
      call c_f_pointer(distances_cptr, distances_fptr, (/k/))
      distances(:) = distances_fptr(:)
    end if
    if (present(lonlats)) then
      call c_f_pointer(lons_cptr, lons_fptr, (/k/))
      call c_f_pointer(lats_cptr, lats_fptr, (/k/))
      lonlats(:,1) = lons_fptr(:)
      lonlats(:,2) = lats_fptr(:)
    end if
    call c_ptr_free(lons_cptr)
    call c_ptr_free(lats_cptr)
    call c_ptr_free(indices_cptr)
    call c_ptr_free(distances_cptr)
  end if
end subroutine IndexKDTree__closestPoints_vectorized_coords

subroutine IndexKDTree__closestPoint_separate_coords(this, plon, plat, index, distance, lon, lat)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  integer(c_int), intent(out) :: index
  real(c_double), intent(out), optional :: distance
  real(c_double), intent(out), optional :: lon
  real(c_double), intent(out), optional :: lat
  real(c_double) :: distance_tmp, lon_tmp, lat_tmp
  call atlas__IndexKDTree__closestPoint(this%CPTR_PGIBUG_A, plon, plat, lon_tmp, lat_tmp, index, distance_tmp)
  if (present(distance)) distance = distance_tmp
  if (present(lon)) lon = lon_tmp
  if (present(lat)) lat = lat_tmp
end subroutine IndexKDTree__closestPoint_separate_coords

subroutine IndexKDTree__closestPoint_vectorized_coords(this, point, index, distance, lonlat)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: point(2)
  integer(c_int), intent(out) :: index
  real(c_double), intent(out), optional :: distance
  real(c_double), intent(out), optional :: lonlat(2)
  real(c_double) :: distance_tmp, lon_tmp, lat_tmp
  call atlas__IndexKDTree__closestPoint(this%CPTR_PGIBUG_A, point(1), point(2), lon_tmp, lat_tmp, index, distance_tmp)
  if (present(distance)) distance = distance_tmp
  if (present(lonlat)) then
     lonlat(1) = lon_tmp
     lonlat(2) = lat_tmp
  end if
end subroutine IndexKDTree__closestPoint_vectorized_coords

subroutine IndexKDTree__closestPointsWithinRadius_separate_coords(this, plon, plat, radius, k, indices, distances, lons, lats)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  real(c_double), intent(in) :: radius
  integer(c_int), intent(out) :: k
  integer(c_int), allocatable, intent(inout), optional :: indices(:)
  real(c_double), allocatable, intent(inout), optional :: distances(:)
  real(c_double), allocatable, intent(inout), optional :: lons(:)
  real(c_double), allocatable, intent(inout), optional :: lats(:)
  integer(c_size_t) :: k_tmp
  type(c_ptr) :: lons_cptr
  type(c_ptr) :: lats_cptr
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  real(c_double), pointer :: lons_fptr(:)
  real(c_double), pointer :: lats_fptr(:)
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  call atlas__IndexKDTree__closestPointsWithinRadius(this%CPTR_PGIBUG_A, plon, plat, radius, &
                                                   & k_tmp, lons_cptr, lats_cptr, indices_cptr, distances_cptr)
  k = int(k_tmp, c_int)
  if (present(indices)) then
    if (allocated(indices)) deallocate(indices)
    allocate(indices(k))
  end if
  if (present(distances)) then
    if (allocated(distances)) deallocate(distances)
    allocate(distances(k))
  end if
  if (present(lons)) then
    if (allocated(lons)) deallocate(lons)
    allocate(lons(k))
  end if
  if (present(lats)) then
    if (allocated(lats)) deallocate(lats)
    allocate(lats(k))
  end if
  if ( k > 0 ) then
    if (present(indices)) then
      call c_f_pointer(indices_cptr, indices_fptr, (/k/))
      indices(:) = indices_fptr(:)
    end if
    if (present(distances)) then
      call c_f_pointer(distances_cptr, distances_fptr, (/k/))
      distances(:) = distances_fptr(:)
    end if
    if (present(lons)) then
      call c_f_pointer(lons_cptr, lons_fptr, (/k/))
      lons(:) = lons_fptr(:)
    end if
    if (present(lats)) then
      call c_f_pointer(lats_cptr, lats_fptr, (/k/))
      lats(:) = lats_fptr(:)
    end if
  end if
  call c_ptr_free(lons_cptr)
  call c_ptr_free(lats_cptr)
  call c_ptr_free(indices_cptr)
  call c_ptr_free(distances_cptr)
end subroutine IndexKDTree__closestPointsWithinRadius_separate_coords

subroutine IndexKDTree__closestPointsWithinRadius_vectorized_coords(this, point, radius, k, indices, distances, lonlats)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: point(2)
  real(c_double), intent(in) :: radius
  integer(c_int), intent(out) :: k
  integer(c_int), allocatable, intent(inout), optional :: indices(:)
  real(c_double), allocatable, intent(inout), optional :: distances(:)
  real(c_double), allocatable, intent(inout), optional :: lonlats(:,:)
  integer(c_size_t) :: k_tmp
  type(c_ptr) :: lons_cptr
  type(c_ptr) :: lats_cptr
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  real(c_double), pointer :: lons_fptr(:)
  real(c_double), pointer :: lats_fptr(:)
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  call atlas__IndexKDTree__closestPointsWithinRadius(this%CPTR_PGIBUG_A, point(1), point(2), radius, &
                                                   & k_tmp, lons_cptr, lats_cptr, indices_cptr, distances_cptr)
  k = int(k_tmp, c_int)
  if (present(indices)) then
    if (allocated(indices)) deallocate(indices)
    allocate(indices(k))
  end if
  if (present(distances)) then
    if (allocated(distances)) deallocate(distances)
    allocate(distances(k))
  end if
  if (present(lonlats)) then
    if (allocated(lonlats)) deallocate(lonlats)
    allocate(lonlats(k,2))
  end if
  if ( k > 0 ) then
    if (present(indices)) then
      call c_f_pointer(indices_cptr, indices_fptr, (/k/))
      indices(:) = indices_fptr(:)
    end if
    if (present(distances)) then
      call c_f_pointer(distances_cptr, distances_fptr, (/k/))
      distances(:) = distances_fptr(:)
    end if
    if (present(lonlats)) then
      call c_f_pointer(lons_cptr, lons_fptr, (/k/))
      call c_f_pointer(lats_cptr, lats_fptr, (/k/))
      lonlats(:,1) = lons_fptr(:)
      lonlats(:,2) = lats_fptr(:)
    end if
  end if
  call c_ptr_free(lons_cptr)
  call c_ptr_free(lats_cptr)
  call c_ptr_free(indices_cptr)
  call c_ptr_free(distances_cptr)
end subroutine IndexKDTree__closestPointsWithinRadius_vectorized_coords

function IndexKDTree__geometry(this) result(geometry)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  type(atlas_Geometry) :: geometry
  call geometry%reset_c_ptr( atlas__IndexKDTree__geometry( this%CPTR_PGIBUG_A ) )
  call geometry%return()
end function IndexKDTree__geometry

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
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
#endif

end module atlas_KDTree_module
