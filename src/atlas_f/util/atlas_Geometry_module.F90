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
use fckit_owned_object_module, only : fckit_owned_object
use atlas_point_module

implicit none

private :: fckit_owned_object

public :: atlas_Geometry

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Geometry

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
  procedure :: xyz2lonlat => Geometry__xyz2lonlat
  procedure :: lonlat2xyz => Geometry__lonlat2xyz
  procedure :: distance_lonlat => Geometry__distance_lonlat
  procedure :: distance_xyz => Geometry__distance_xyz
  generic :: distance => distance_lonlat, distance_xyz
  procedure :: radius => Geometry__radius
  procedure :: area => Geometry__area
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Geometry__final_auto
#endif

END TYPE atlas_Geometry

!------------------------------------------------------------------------------

interface atlas_Geometry
  module procedure atlas_Geometry__ctor_name
  module procedure atlas_Geometry__ctor_radius
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! Geometry routines

function atlas_Geometry__ctor_name(name) result(this)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  character(kind=c_char,len=*), intent(in) :: name
  type(atlas_Geometry) :: this
  call this%reset_c_ptr( atlas__Geometry__new_name( c_str(name) ) )
  call this%return()
end function atlas_Geometry__ctor_name

function atlas_Geometry__ctor_radius(radius) result(this)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  real(c_double), intent(in) :: radius
  type(atlas_Geometry) :: this
  call this%reset_c_ptr( atlas__Geometry__new_radius( radius ) )
  call this%return()
end function atlas_Geometry__ctor_radius

subroutine atlas_Geometry__delete(this)
  use atlas_Geometry_c_binding
  class(atlas_Geometry), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Geometry__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_Geometry__delete

subroutine Geometry__xyz2lonlat(this, x, y, z, lon, lat)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: y
  real(c_double), intent(in) :: z
  real(c_double), intent(out) :: lon
  real(c_double), intent(out) :: lat
  call atlas__Geometry__xyz2lonlat(this%CPTR_PGIBUG_A, x, y, z, lon, lat)
end subroutine Geometry__xyz2lonlat

subroutine Geometry__lonlat2xyz(this, lon, lat, x, y, z)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  real(c_double), intent(out) :: x
  real(c_double), intent(out) :: y
  real(c_double), intent(out) :: z
  call atlas__Geometry__lonlat2xyz(this%CPTR_PGIBUG_A, lon, lat, x, y, z)
end subroutine Geometry__lonlat2xyz

function Geometry__distance_lonlat(this, lon1, lat1, lon2, lat2) result(distance_lonlat)
  use atlas_Geometry_c_binding
  use fckit_c_interop_module
  class(atlas_Geometry), intent(in) :: this
  real(c_double), intent(in) :: lon1
  real(c_double), intent(in) :: lat1
  real(c_double), intent(in) :: lon2
  real(c_double), intent(in) :: lat2
  real(c_double) :: distance_lonlat
  distance_lonlat = atlas__Geometry__distance_lonlat(this%CPTR_PGIBUG_A, lon1, lat1, lon2, lat2)
end function Geometry__distance_lonlat

function Geometry__distance_xyz(this, x1, y1, z1, x2, y2, z2) result(distance_xyz)
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
  distance_xyz = atlas__Geometry__distance_xyz(this%CPTR_PGIBUG_A, x1, y1, z1, x2, y2, z2)
end function Geometry__distance_xyz

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
