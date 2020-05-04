! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_point_module

use, intrinsic :: iso_c_binding
use fckit_object_module, only : fckit_object

implicit none

private :: fckit_object

public :: atlas_PointXY, atlas_PointXYZ, atlas_PointLonLat

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_PointXY

! Purpose :
! -------
!   *PointXY* : Container of PointXY

! Methods :
! -------

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_PointXY__delete
  procedure :: print => PointXY__print
  procedure :: assign => PointXY__assign
  procedure :: x => PointXY__x
  procedure :: y => PointXY__y
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_PointXY__final_auto
#endif

END TYPE atlas_PointXY

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_PointXYZ

! Purpose :
! -------
!   *PointXYZ* : Container of PointXYZ

! Methods :
! -------

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_PointXYZ__delete
  procedure :: print => PointXYZ__print
  procedure :: assign => PointXYZ__assign
  procedure :: x => PointXYZ__x
  procedure :: y => PointXYZ__y
  procedure :: z => PointXYZ__z
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_PointXYZ__final_auto
#endif

END TYPE atlas_PointXYZ

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_PointLonLat

! Purpose :
! -------
!   *PointLonLat* : Container of PointLonLat

! Methods :
! -------
!   assign : assign coordinates
!   print : print coordinates

! Author :
! ------
!   April-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_PointLonLat__delete
  procedure :: print => PointLonLat__print
  procedure :: assign => PointLonLat__assign
  procedure :: lon => PointLonLat__lon
  procedure :: lat => PointLonLat__lat
  procedure :: normalise_only => PointLonLat__normalise_only
  procedure :: normalise_west => PointLonLat__normalise_west
  procedure :: normalise_west_east => PointLonLat__normalise_west_east
  generic :: normalise => normalise_only, normalise_west, normalise_west_east
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_PointLonLat__final_auto
#endif

END TYPE atlas_PointLonLat

!------------------------------------------------------------------------------

interface atlas_PointXY
  module procedure atlas_PointXY__cptr
  module procedure atlas_PointXY__ctor
  module procedure atlas_PointXY__ctor_xy
end interface

interface atlas_PointXYZ
  module procedure atlas_PointXYZ__cptr
  module procedure atlas_PointXYZ__ctor
  module procedure atlas_PointXYZ__ctor_xyz
end interface

interface atlas_PointLonLat
  module procedure atlas_PointLonLat__cptr
  module procedure atlas_PointLonLat__ctor
  module procedure atlas_PointLonLat__ctor_lonlat
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! PointXY routines

function atlas_PointXY__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_PointXY) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function atlas_PointXY__cptr

function atlas_PointXY__ctor() result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointXY) :: this
  call this%reset_c_ptr( atlas__PointXY__new() )
end function atlas_PointXY__ctor

subroutine atlas_PointXY__delete(this)
  use atlas_point_c_binding
  class(atlas_PointXY), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__PointXY__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_PointXY__delete

subroutine PointXY__print(this, channel)
  use atlas_point_c_binding
  use fckit_log_module, only : fckit_logchannel
  class(atlas_PointXY), intent(in) :: this
  type(fckit_logchannel), intent(in) :: channel
  call atlas__PointXY__print(this%CPTR_PGIBUG_A,channel%CPTR_PGIBUG_A)
end subroutine PointXY__print

function atlas_PointXY__ctor_xy(x, y) result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointXY) :: this
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: y
  call this%reset_c_ptr( atlas__PointXY__new_xy(x, y) )
end function atlas_PointXY__ctor_xy

function PointXY__x(this) result(x)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointXY), intent(in) :: this
  real(c_double) :: x
  x = atlas__PointXY__x(this%CPTR_PGIBUG_A)
end function PointXY__x

function PointXY__y(this) result(y)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointXY), intent(in) :: this
  real(c_double) :: y
  y = atlas__PointXY__y(this%CPTR_PGIBUG_A)
end function PointXY__y

subroutine PointXY__assign(this, x, y)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointXY), intent(inout) :: this
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: y
  call atlas__PointXY__assign(this%CPTR_PGIBUG_A, x, y)
end subroutine PointXY__assign

! -----------------------------------------------------------------------------
! PointXYZ routines

function atlas_PointXYZ__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_PointXYZ) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function atlas_PointXYZ__cptr

function atlas_PointXYZ__ctor() result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointXYZ) :: this
  call this%reset_c_ptr( atlas__PointXYZ__new() )
end function atlas_PointXYZ__ctor

subroutine atlas_PointXYZ__delete(this)
  use atlas_point_c_binding
  class(atlas_PointXYZ), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__PointXYZ__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_PointXYZ__delete

subroutine PointXYZ__print(this, channel)
  use atlas_point_c_binding
  use fckit_log_module, only : fckit_logchannel
  class(atlas_PointXYZ), intent(in) :: this
  type(fckit_logchannel), intent(in) :: channel
  call atlas__PointXYZ__print(this%CPTR_PGIBUG_A,channel%CPTR_PGIBUG_A)
end subroutine PointXYZ__print

function atlas_PointXYZ__ctor_XYZ(x, y, z) result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointXYZ) :: this
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: y
  real(c_double), intent(in) :: z
  call this%reset_c_ptr( atlas__PointXYZ__new_XYZ(x, y, z) )
end function atlas_PointXYZ__ctor_XYZ

function PointXYZ__x(this) result(x)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointXYZ), intent(in) :: this
  real(c_double) :: x
  x = atlas__PointXYZ__x(this%CPTR_PGIBUG_A)
end function PointXYZ__x

function PointXYZ__y(this) result(y)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointXYZ), intent(in) :: this
  real(c_double) :: y
  y = atlas__PointXYZ__y(this%CPTR_PGIBUG_A)
end function PointXYZ__y

function PointXYZ__z(this) result(z)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointXYZ), intent(in) :: this
  real(c_double) :: z
  z = atlas__PointXYZ__z(this%CPTR_PGIBUG_A)
end function PointXYZ__z

subroutine PointXYZ__assign(this, x, y, z)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointXYZ), intent(inout) :: this
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: y
  real(c_double), intent(in) :: z
  call atlas__PointXYZ__assign(this%CPTR_PGIBUG_A, x, y, z)
end subroutine PointXYZ__assign

! -----------------------------------------------------------------------------
! PointLonLat routines

function atlas_PointLonLat__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_PointLonLat) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function atlas_PointLonLat__cptr

function atlas_PointLonLat__ctor() result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointLonLat) :: this
  call this%reset_c_ptr( atlas__PointLonLat__new() )
end function atlas_PointLonLat__ctor

subroutine atlas_PointLonLat__delete(this)
  use atlas_point_c_binding
  class(atlas_PointLonLat), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__PointLonLat__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_PointLonLat__delete

subroutine PointLonLat__print(this, channel)
  use atlas_point_c_binding
  use fckit_log_module, only : fckit_logchannel
  class(atlas_PointLonLat), intent(in) :: this
  type(fckit_logchannel), intent(in) :: channel
  call atlas__PointLonLat__print(this%CPTR_PGIBUG_A,channel%CPTR_PGIBUG_A)
end subroutine PointLonLat__print

function atlas_PointLonLat__ctor_lonlat(lon, lat) result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointLonLat) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  call this%reset_c_ptr( atlas__PointLonLat__new_lonlat(lon, lat) )
end function atlas_PointLonLat__ctor_lonlat

function PointLonLat__lon(this) result(lon)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointLonLat), intent(in) :: this
  real(c_double) :: lon
  lon = atlas__PointLonLat__lon(this%CPTR_PGIBUG_A)
end function PointLonLat__lon

function PointLonLat__lat(this) result(lat)
  use atlas_point_c_binding
  use fckit_c_interop_module
  class(atlas_PointLonLat), intent(in) :: this
  real(c_double) :: lat
  lat = atlas__PointLonLat__lat(this%CPTR_PGIBUG_A)
end function PointLonLat__lat

subroutine PointLonLat__assign(this, lon, lat)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointLonLat), intent(inout) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  call atlas__PointLonLat__assign(this%CPTR_PGIBUG_A, lon, lat)
end subroutine PointLonLat__assign

subroutine PointLonLat__normalise_only(this)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointLonLat), intent(inout) :: this
  call atlas__PointLonLat__normalise(this%CPTR_PGIBUG_A)
end subroutine PointLonLat__normalise_only

subroutine PointLonLat__normalise_west(this, west)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointLonLat), intent(inout) :: this
  real(c_double), intent(in) :: west
  call atlas__PointLonLat__normalise_west(this%CPTR_PGIBUG_A, west)
end subroutine PointLonLat__normalise_west

subroutine PointLonLat__normalise_west_east(this, west, east)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointLonLat), intent(inout) :: this
  real(c_double), intent(in) :: west
  real(c_double), intent(in) :: east
  call atlas__PointLonLat__normalise_west_east(this%CPTR_PGIBUG_A, west, east)
end subroutine PointLonLat__normalise_west_east

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_PointXY__final_auto(this)
  type(atlas_PointXY), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_PointXY__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_PointXYZ__final_auto(this)
  type(atlas_PointXYZ), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_PointXYZ__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_PointLonLat__final_auto(this)
  type(atlas_PointLonLat), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_PointLonLat__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

end module atlas_point_module
