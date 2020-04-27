! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_earth_module

use, intrinsic :: iso_c_binding
use fckit_object_module, only : fckit_object
use atlas_point_module

implicit none

private :: fckit_object

public :: atlas_Earth

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_Earth

! Purpose :
! -------
!   *Earth* : Container of Earth

! Methods :
! -------

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_Earth__delete
  procedure :: radius => Earth__radius
  procedure :: central_angle_lonlat => Earth__central_angle_lonlat
  procedure :: central_angle_xyz => Earth__central_angle_xyz
  generic :: central_angle => central_angle_lonlat, central_angle_xyz
  procedure :: distance_lonlat => Earth__distance_lonlat
  procedure :: distance_xyz => Earth__distance_xyz
  generic :: distance => distance_lonlat, distance_xyz
  procedure :: area_total => Earth__area
  procedure :: area_wn_es => Earth__area_wn_es
  generic :: area => area_total, area_wn_es
  procedure :: great_circle_latitude_given_longitude => Earth__great_circle_latitude_given_longitude
  procedure :: great_circle_longitude_given_latitude => Earth__great_circle_longitude_given_latitude
  procedure :: convert_spherical_to_cartesian => Earth__convert_spherical_to_cartesian
  procedure :: convert_cartesian_to_spherical => Earth__convert_cartesian_to_spherical
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Earth__final_auto
#endif

END TYPE atlas_Earth

!------------------------------------------------------------------------------

interface atlas_Earth
  module procedure atlas_Earth__ctor
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! Earth routines

function atlas_Earth__ctor() result(this)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  type(atlas_Earth) :: this
  call this%reset_c_ptr( atlas__Earth__new() )
end function atlas_Earth__ctor

subroutine atlas_Earth__delete(this)
  use atlas_earth_c_binding
  class(atlas_Earth), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Earth__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_Earth__delete

function Earth__radius(this) result(radius)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  real(c_double) :: radius
  radius = atlas__Earth__radius(this%CPTR_PGIBUG_A)
end function Earth__radius

function Earth__central_angle_lonlat(this, Alonlat, Blonlat) result(central_angle_lonlat)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double) :: central_angle_lonlat
  central_angle_lonlat = atlas__Earth__central_angle_2(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A)
end function Earth__central_angle_lonlat

function Earth__central_angle_xyz(this, A, B) result(central_angle_xyz)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: A
  type(atlas_PointXYZ), intent(in) :: B
  real(c_double) :: central_angle_xyz
  central_angle_xyz = atlas__Earth__central_angle_3(this%CPTR_PGIBUG_A, A%CPTR_PGIBUG_A, B%CPTR_PGIBUG_A)
end function Earth__central_angle_xyz

function Earth__distance_lonlat(this, Alonlat, Blonlat) result(distance_lonlat)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double) :: distance_lonlat
  distance_lonlat = atlas__Earth__distance_2(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A)
end function Earth__distance_lonlat

function Earth__distance_xyz(this, A, B) result(distance_xyz)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: A
  type(atlas_PointXYZ), intent(in) :: B
  real(c_double) :: distance_xyz
  distance_xyz = atlas__Earth__distance_3(this%CPTR_PGIBUG_A, A%CPTR_PGIBUG_A, B%CPTR_PGIBUG_A)
end function Earth__distance_xyz

function Earth__area(this) result(area)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  real(c_double) :: area
  area = atlas__Earth__area(this%CPTR_PGIBUG_A)
end function Earth__area

function Earth__area_wn_es(this, WestNorth, EastSouth) result(area_wn_es)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: WestNorth
  type(atlas_PointLonLat), intent(in) :: EastSouth
  real(c_double) :: area_wn_es
  area_wn_es = atlas__Earth__area_wn_es(this%CPTR_PGIBUG_A, WestNorth%CPTR_PGIBUG_A, EastSouth%CPTR_PGIBUG_A)
end function Earth__area_wn_es

function Earth__great_circle_latitude_given_longitude(this, Alonlat, Blonlat, Clon) result(great_circle_latitude_given_longitude)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double), intent(in) :: Clon
  real(c_double) :: great_circle_latitude_given_longitude
  great_circle_latitude_given_longitude = atlas__Earth__great_circle_latitude_given_longitude(this%CPTR_PGIBUG_A, &
 & Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A, Clon)
end function Earth__great_circle_latitude_given_longitude

subroutine Earth__great_circle_longitude_given_latitude(this, Alonlat, Blonlat, Clat, Clon1, Clon2)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointLonLat), intent(in) :: Blonlat
  real(c_double), intent(in) :: Clat
  real(c_double), intent(inout) :: Clon1
  real(c_double), intent(inout) :: Clon2
  call atlas__Earth__great_circle_longitude_given_latitude(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A, &
 & Clat, Clon1, Clon2)
end subroutine Earth__great_circle_longitude_given_latitude

subroutine Earth__convert_spherical_to_cartesian(this, Alonlat, B)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: Alonlat
  type(atlas_PointXYZ), intent(inout) :: B
  call atlas__Earth__convert_spherical_to_cartesian(this%CPTR_PGIBUG_A, Alonlat%CPTR_PGIBUG_A, B%CPTR_PGIBUG_A)
end subroutine Earth__convert_spherical_to_cartesian

subroutine Earth__convert_cartesian_to_spherical(this, A, Blonlat)
  use atlas_earth_c_binding
  use fckit_c_interop_module
  class(atlas_Earth), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: A
  type(atlas_PointLonLat), intent(inout) :: Blonlat
  call atlas__Earth__convert_cartesian_to_spherical(this%CPTR_PGIBUG_A, A%CPTR_PGIBUG_A, Blonlat%CPTR_PGIBUG_A)
end subroutine Earth__convert_cartesian_to_spherical

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_Earth__final_auto(this)
  type(atlas_Earth), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Earth__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

end module atlas_earth_module
