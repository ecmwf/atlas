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

public :: atlas_PointLonLat

private

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
  procedure :: normalise_ => PointLonLat__normalise
  procedure :: normalise_west => PointLonLat__normalise_west
  procedure :: normalise_west_east => PointLonLat__normalise_west_east
  generic :: normalise => normalise_, normalise_west, normalise_west_east
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_PointLonLat__final_auto
#endif

END TYPE atlas_PointLonLat

!------------------------------------------------------------------------------

interface atlas_PointLonLat
  module procedure atlas_PointLonLat__ctor
  module procedure atlas_PointLonLat__ctor_lonlat
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! PointLonLat routines

function atlas_PointLonLat__ctor() result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointLonLat) :: this
  call this%reset_c_ptr( atlas__PointLonLat__new() )
end function atlas_PointLonLat__ctor

function atlas_PointLonLat__ctor_lonlat(lon, lat) result(this)
  use atlas_point_c_binding
  use fckit_c_interop_module
  type(atlas_PointLonLat) :: this
  real(c_double), intent(in) :: lon
  real(c_double), intent(in) :: lat
  call this%reset_c_ptr( atlas__PointLonLat__new_lonlat(lon, lat) )
end function atlas_PointLonLat__ctor_lonlat

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

subroutine PointLonLat__normalise(this)
  use atlas_point_c_binding
  use, intrinsic :: iso_c_binding
  class(atlas_PointLonLat), intent(inout) :: this
  call atlas__PointLonLat__normalise(this%CPTR_PGIBUG_A)
end subroutine PointLonLat__normalise

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
