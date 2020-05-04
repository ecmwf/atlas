! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Geometry_module

use, intrinsic :: iso_c_binding
use fckit_object_module, only : fckit_object
use atlas_point_module

implicit none

private :: fckit_object

public :: atlas_Geometry

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_Geometry

! Purpose :
! -------
!   *Geometry* : Container of Geometry

! Methods :
! -------

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_Geometry__delete
  procedure :: xyz2lonlat_point => Geometry__xyz2lonlat_point
  procedure :: xyz2lonlat_real => Geometry__xyz2lonlat_real
  generic :: xyz2lonlat => xyz2lonlat_point, xyz2lonlat_real
  procedure :: lonlat2xyz_point => Geometry__lonlat2xyz_point
  procedure :: lonlat2xyz_real => Geometry__lonlat2xyz_real
  generic :: lonlat2xyz => lonlat2xyz_point, lonlat2xyz_real
  procedure :: distance_lonlat_point => Geometry__distance_lonlat_point
  procedure :: distance_lonlat_real => Geometry__distance_lonlat_real
  procedure :: distance_xyz_point => Geometry__distance_xyz_point
  procedure :: distance_xyz_real => Geometry__distance_xyz_real
  generic :: distance => distance_lonlat_point, distance_lonlat_real, &
                       & distance_xyz_point, distance_xyz_real
  procedure :: radius => Geometry__radius
  procedure :: area => Geometry__area
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Geometry__final_auto
#endif

END TYPE atlas_Geometry

!------------------------------------------------------------------------------

interface atlas_Geometry
  module procedure atlas_Geometry__cptr
  module procedure atlas_Geometry__ctor_name
  module procedure atlas_Geometry__ctor_radius
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! Geometry routines

function atlas_Geometry__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_Geometry) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function atlas_Geometry__cptr

function atlas_Geometry__ctor_name(name) result(this)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  character(kind=c_char,len=*), intent(in) :: name
  type(atlas_Geometry) :: this
  call this%reset_c_ptr( atlas__Geometry__new_name( c_str(name) ) )
end function atlas_Geometry__ctor_name

function atlas_Geometry__ctor_radius(radius) result(this)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  real(c_double), intent(in) :: radius
  type(atlas_Geometry) :: this
  call this%reset_c_ptr( atlas__Geometry__new_radius( radius ) )
end function atlas_Geometry__ctor_radius

subroutine atlas_Geometry__delete(this)
  use atlas_Geometry_c_binding
  class(atlas_Geometry), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Geometry__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_Geometry__delete

subroutine Geometry__xyz2lonlat_point(this, xyz, lonlat)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: xyz
  type(atlas_PointLonLat), intent(inout) :: lonlat
  call atlas__Geometry__xyz2lonlat(this%CPTR_PGIBUG_A, xyz%CPTR_PGIBUG_A, lonlat%CPTR_PGIBUG_A)
end subroutine Geometry__xyz2lonlat_point

subroutine Geometry__xyz2lonlat_real(this, x, y, z, lon, lat)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: y
  real(c_double), intent(in) :: z
  real(c_double), intent(out) :: lon
  real(c_double), intent(out) :: lat
  type(atlas_PointXYZ) :: xyz
  type(atlas_PointLonLat) :: lonlat
  xyz = atlas_PointXYZ(x, y, z)
  lonlat = atlas_PointLonLat()
  call this%xyz2lonlat(xyz, lonlat)
  lon = lonlat%lon()
  lat = lonlat%lat()
end subroutine Geometry__xyz2lonlat_real

subroutine Geometry__lonlat2xyz_point(this, lonlat, xyz)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: lonlat
  type(atlas_PointXYZ), intent(inout) :: xyz
  type(atlas_PointLonLat) :: lonlat_tmp
  lonlat_tmp = atlas_PointLonLat(lonlat%lon(), lonlat%lat())
  call lonlat_tmp%normalise()
  call atlas__Geometry__lonlat2xyz(this%CPTR_PGIBUG_A, lonlat_tmp%CPTR_PGIBUG_A, xyz%CPTR_PGIBUG_A)
end subroutine Geometry__lonlat2xyz_point

subroutine Geometry__lonlat2xyz_real(this, lon, lat, x, y, z)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  real(c_double), intent(out) :: x
  real(c_double), intent(out) :: y
  real(c_double), intent(out) :: z
  type(atlas_PointLonLat) :: lonlat
  type(atlas_PointXYZ) :: xyz
  lonlat = atlas_PointLonLat(lon, lat)
  xyz = atlas_PointXYZ()
  call this%lonlat2xyz(lonlat, xyz)
  x = xyz%x()
  y = xyz%y()
  z = xyz%z()
end subroutine Geometry__lonlat2xyz_real

function Geometry__distance_lonlat_point(this, p1, p2) result(distance_lonlat)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  type(atlas_PointLonLat), intent(in) :: p1
  type(atlas_PointLonLat), intent(in) :: p2
  real(c_double) :: distance_lonlat
  type(atlas_PointLonLat) :: p1_tmp, p2_tmp
  p1_tmp = atlas_PointLonLat(p1%lon(), p1%lat())
  p2_tmp = atlas_PointLonLat(p2%lon(), p2%lat())
  call p1_tmp%normalise()
  call p2_tmp%normalise()
  distance_lonlat = atlas__Geometry__distance_2(this%CPTR_PGIBUG_A, p1_tmp%CPTR_PGIBUG_A, p2_tmp%CPTR_PGIBUG_A)
end function Geometry__distance_lonlat_point

function Geometry__distance_lonlat_real(this, lon1, lat1, lon2, lat2) result(distance_lonlat)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: lon1
  real(c_double), intent(in) :: lat1
  real(c_double), intent(in) :: lon2
  real(c_double), intent(in) :: lat2
  real(c_double) :: distance_lonlat
  type(atlas_PointLonLat) :: p1
  type(atlas_PointLonLat) :: p2
  p1 = atlas_PointLonLat(lon1, lat1)
  p2 = atlas_PointLonLat(lon2, lat2)
  distance_lonlat = this%distance(p1, p2)
end function Geometry__distance_lonlat_real

function Geometry__distance_xyz_point(this, p1, p2) result(distance_xyz)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  type(atlas_PointXYZ), intent(in) :: p1
  type(atlas_PointXYZ), intent(in) :: p2
  real(c_double) :: distance_xyz
  distance_xyz = atlas__Geometry__distance_3(this%CPTR_PGIBUG_A, p1%CPTR_PGIBUG_A, p2%CPTR_PGIBUG_A)
end function Geometry__distance_xyz_point

function Geometry__distance_xyz_real(this, x1, y1, z1, x2, y2, z2) result(distance_xyz)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: x1
  real(c_double), intent(in) :: y1
  real(c_double), intent(in) :: z1
  real(c_double), intent(in) :: x2
  real(c_double), intent(in) :: y2
  real(c_double), intent(in) :: z2
  real(c_double) :: distance_xyz
  type(atlas_PointXYZ) :: p1
  type(atlas_PointXYZ) :: p2
  p1 = atlas_PointXYZ(x1, y1, z1)
  p2 = atlas_PointXYZ(x2, y2, z2)
  distance_xyz = this%distance(p1, p2)
end function Geometry__distance_xyz_real

function Geometry__radius(this) result(radius)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double) :: radius
  radius = atlas__Geometry__radius(this%CPTR_PGIBUG_A)
end function Geometry__radius

function Geometry__area(this) result(area)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double) :: area
  area = atlas__Geometry__area(this%CPTR_PGIBUG_A)
end function Geometry__area

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_Geometry__final_auto(this)
  type(atlas_Geometry), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Geometry__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

end module atlas_Geometry_module
