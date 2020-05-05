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
  procedure :: insert => IndexKDTree__insert
  procedure :: build_only => IndexKDTree__build_only
  procedure :: build_list => IndexKDTree__build_list
  generic :: build => build_only, build_list
  procedure :: closestPoints => IndexKDTree__closestPoints
  procedure :: closestPoint => IndexKDTree__closestPoint
  procedure :: closestPointsWithinRadius => IndexKDTree__closestPointsWithinRadius
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

subroutine IndexKDTree__insert(this, lon, lat, index)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  integer(c_int), intent(in) :: index
  call atlas__IndexKDTree__insert(this%CPTR_PGIBUG_A, lon, lat, index)
end subroutine IndexKDTree__insert

subroutine IndexKDTree__build_only(this)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  call atlas__IndexKDTree__build(this%CPTR_PGIBUG_A)
end subroutine IndexKDTree__build_only

subroutine IndexKDTree__build_list(this, k, lons, lats, indices)
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
end subroutine IndexKDTree__build_list

subroutine IndexKDTree__closestPoints(this, plon, plat, k, lons, lats, indices, distances)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  integer(c_int), intent(in) :: k
  real(c_double), intent(out) :: lons(k)
  real(c_double), intent(out) :: lats(k)
  integer(c_int), intent(out) :: indices(k)
  real(c_double), intent(out) :: distances(k)
  type(c_ptr) :: lons_cptr
  type(c_ptr) :: lats_cptr
  type(c_ptr) :: indices_cptr
  type(c_ptr) :: distances_cptr
  real(c_double), pointer :: lons_fptr(:)
  real(c_double), pointer :: lats_fptr(:)
  integer(c_int), pointer :: indices_fptr(:)
  real(c_double), pointer :: distances_fptr(:)
  if ( k > 0 ) then
    call atlas__IndexKDTree__closestPoints(this%CPTR_PGIBUG_A, plon, plat, int(k, c_size_t), lons_cptr, lats_cptr, indices_cptr, &
                                       & distances_cptr)
    call c_f_pointer(lons_cptr, lons_fptr, (/k/))
    call c_f_pointer(lats_cptr, lats_fptr, (/k/))
    call c_f_pointer(indices_cptr, indices_fptr, (/k/))
    call c_f_pointer(distances_cptr, distances_fptr, (/k/))
    lons(:) = lons_fptr(:)
    lats(:) = lats_fptr(:)
    indices(:) = indices_fptr(:)
    distances(:) = distances_fptr(:)
    call c_ptr_free(lons_cptr)
    call c_ptr_free(lats_cptr)
    call c_ptr_free(indices_cptr)
    call c_ptr_free(distances_cptr)
  end if
end subroutine IndexKDTree__closestPoints

subroutine IndexKDTree__closestPoint(this, plon, plat, lon, lat, index, distance)
  use atlas_KDTree_c_binding
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  real(c_double), intent(out) :: lon
  real(c_double), intent(out) :: lat
  integer(c_int), intent(out) :: index
  real(c_double), intent(out) :: distance
  call atlas__IndexKDTree__closestPoint(this%CPTR_PGIBUG_A, plon, plat, lon, lat, index, distance)
end subroutine IndexKDTree__closestPoint

subroutine IndexKDTree__closestPointsWithinRadius(this, plon, plat, radius, k, lons, lats, indices, distances)
  use atlas_KDTree_c_binding
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  use fckit_c_interop_module
  class(atlas_IndexKDTree), intent(in) :: this
  real(c_double), intent(in) :: plon
  real(c_double), intent(in) :: plat
  real(c_double), intent(in) :: radius
  integer(c_int), intent(out) :: k
  real(c_double), allocatable, intent(inout) :: lons(:)
  real(c_double), allocatable, intent(inout) :: lats(:)
  integer(c_int), allocatable, intent(inout) :: indices(:)
  real(c_double), allocatable, intent(inout) :: distances(:)
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
  if (allocated(lons)) deallocate(lons)
  if (allocated(lats)) deallocate(lats)
  if (allocated(indices)) deallocate(indices)
  if (allocated(distances)) deallocate(distances)
  allocate(lons(k))
  allocate(lats(k))
  allocate(indices(k))
  allocate(distances(k))
  if ( k > 0 ) then
    call c_f_pointer(lons_cptr, lons_fptr, (/k/))
    call c_f_pointer(lats_cptr, lats_fptr, (/k/))
    call c_f_pointer(indices_cptr, indices_fptr, (/k/))
    call c_f_pointer(distances_cptr, distances_fptr, (/k/))
    lons(:) = lons_fptr(:)
    lats(:) = lats_fptr(:)
    indices(:) = indices_fptr(:)
    distances(:) = distances_fptr(:)
 end if
 call c_ptr_free(lons_cptr)
 call c_ptr_free(lats_cptr)
 call c_ptr_free(indices_cptr)
 call c_ptr_free(distances_cptr)
end subroutine IndexKDTree__closestPointsWithinRadius

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
