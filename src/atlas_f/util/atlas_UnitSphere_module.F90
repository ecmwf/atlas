! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_UnitSphere_module

use, intrinsic :: iso_c_binding
use fckit_object_module, only : fckit_object
use atlas_point_module

implicit none

private :: fckit_object

public :: atlas_UnitSphere

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_UnitSphere

! Purpose :
! -------
!   *UnitSphere* : Container of UnitSphere

! Methods :
! -------

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_UnitSphere__delete
  procedure :: radius => UnitSphere__radius
  procedure :: central_angle_lonlat => UnitSphere__central_angle_lonlat
  procedure :: central_angle_xyz => UnitSphere__central_angle_xyz
  generic :: central_angle => central_angle_lonlat, central_angle_xyz
  procedure :: distance_lonlat => UnitSphere__distance_lonlat
  procedure :: distance_xyz => UnitSphere__distance_xyz
  generic :: distance => distance_lonlat, distance_xyz
  procedure :: area_total => UnitSphere__area
  procedure :: area_wn_es => UnitSphere__area_wn_es
  generic :: area => area_total, area_wn_es
  procedure :: great_circle_latitude_given_longitude => UnitSphere__great_circle_latitude_given_longitude
  procedure :: great_circle_longitude_given_latitude => UnitSphere__great_circle_longitude_given_latitude
  procedure :: convert_spherical_to_cartesian => UnitSphere__convert_spherical_to_cartesian
  procedure :: convert_cartesian_to_spherical => UnitSphere__convert_cartesian_to_spherical
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_UnitSphere__final_auto
#endif

END TYPE atlas_UnitSphere

!------------------------------------------------------------------------------

interface atlas_UnitSphere
  module procedure atlas_UnitSphere__ctor
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! UnitSphere routines

function atlas_UnitSphere__ctor() result(this)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  type(atlas_UnitSphere) :: this
  call this%reset_c_ptr( atlas__UnitSphere__new() )
end function atlas_UnitSphere__ctor

subroutine atlas_UnitSphere__delete(this)
  use atlas_UnitSphere_c_binding
  class(atlas_UnitSphere), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__UnitSphere__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_UnitSphere__delete

function UnitSphere__radius(this) result(radius)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  real(c_double) :: radius
  radius = atlas__UnitSphere__radius(this%CPTR_PGIBUG_A)
end function UnitSphere__radius

function UnitSphere__central_angle_lonlat(this, Alonlat, Blonlat) result(central_angle_lonlat)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double) :: central_angle_lonlat
  central_angle_lonlat = atlas__UnitSphere__central_angle_2(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A)
end function UnitSphere__central_angle_lonlat

function UnitSphere__central_angle_xyz(this, A, B) result(central_angle_xyz)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: A
  type(atlas_PointXYZ), intent(in) :: B
  real(c_double) :: central_angle_xyz
  central_angle_xyz = atlas__UnitSphere__central_angle_3(this%CPTR_PGIBUG_A, A%CPTR_PGIBUG_A, B%CPTR_PGIBUG_A)
end function UnitSphere__central_angle_xyz

function UnitSphere__distance_lonlat(this, Alonlat, Blonlat) result(distance_lonlat)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double) :: distance_lonlat
  distance_lonlat = atlas__UnitSphere__distance_2(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A)
end function UnitSphere__distance_lonlat

function UnitSphere__distance_xyz(this, A, B) result(distance_xyz)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: A
  type(atlas_PointXYZ), intent(in) :: B
  real(c_double) :: distance_xyz
  distance_xyz = atlas__UnitSphere__distance_3(this%CPTR_PGIBUG_A, A%CPTR_PGIBUG_A, B%CPTR_PGIBUG_A)
end function UnitSphere__distance_xyz

function UnitSphere__area(this) result(area)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  real(c_double) :: area
  area = atlas__UnitSphere__area(this%CPTR_PGIBUG_A)
end function UnitSphere__area

function UnitSphere__area_wn_es(this, WestNorth, EastSouth) result(area_wn_es)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: WestNorth
  type(atlas_PointLonLat), intent(in) :: EastSouth
  real(c_double) :: area_wn_es
  area_wn_es = atlas__UnitSphere__area_wn_es(this%CPTR_PGIBUG_A, WestNorth%CPTR_PGIBUG_A, EastSouth%CPTR_PGIBUG_A)
end function UnitSphere__area_wn_es

function UnitSphere__great_circle_latitude_given_longitude(this, Alonlat, Blonlat, Clon) &
 & result(great_circle_latitude_given_longitude)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double), intent(in) :: Clon
  real(c_double) :: great_circle_latitude_given_longitude
  great_circle_latitude_given_longitude = atlas__UnitSphere__great_circle_latitude_given_longitude(this%CPTR_PGIBUG_A, &
 & Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A, Clon)
end function UnitSphere__great_circle_latitude_given_longitude

subroutine UnitSphere__great_circle_longitude_given_latitude(this, Alonlat, Blonlat, Clat, Clon1, Clon2)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double), intent(in) :: Clat
  real(c_double), intent(inout) :: Clon1
  real(c_double), intent(inout) :: Clon2
  call atlas__UnitSphere__great_circle_longitude_given_latitude(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, &
 & Blonlat%CPTR_PGIBUG_A, Clat, Clon1, Clon2)
end subroutine UnitSphere__great_circle_longitude_given_latitude

subroutine UnitSphere__convert_spherical_to_cartesian(this, Alonlat, B)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointXYZ), intent(inout) :: B
  call atlas__UnitSphere__convert_spherical_to_cartesian(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, B%CPTR_PGIBUG_A)
end subroutine UnitSphere__convert_spherical_to_cartesian

subroutine UnitSphere__convert_cartesian_to_spherical(this, A, Blonlat)
  use atlas_UnitSphere_c_binding
  use fckit_c_interop_module
  class(atlas_UnitSphere), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: A
  type(atlas_PointLonLat), intent(inout) :: Blonlat
  call atlas__UnitSphere__convert_cartesian_to_spherical(this%CPTR_PGIBUG_A, A%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A)
end subroutine UnitSphere__convert_cartesian_to_spherical

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_UnitSphere__final_auto(this)
  type(atlas_UnitSphere), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_UnitSphere__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

end module atlas_UnitSphere_module
