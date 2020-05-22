! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Projection_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_config_module, only : atlas_Config

implicit none

private :: fckit_owned_object
private :: atlas_Config

public :: atlas_Projection
public :: atlas_RotatedLonLatProjection
public :: atlas_LambertConformalConicProjection
public :: atlas_RotatedSchmidtProjection

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Projection

! Purpose :
! -------
!   *atlas_Projection* :

! Methods :
! -------
!   type: The name or tag this function space was created with

! Author :
! ------
!   May-2020 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: type => atlas_Projection__type
  procedure, public :: hash
  procedure, public :: spec

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Projection__final_auto
#endif

END TYPE atlas_Projection

interface atlas_Projection
  module procedure atlas_Projection__ctor_cptr
  module procedure atlas_Projection__ctor_config
end interface

interface atlas_RotatedSchmidtProjection
  module procedure atlas_RotatedSchmidtProjection_real64
end interface

interface atlas_LambertConformalConicProjection
  module procedure atlas_LambertConformalConicProjection_real64
end interface

interface atlas_RotatedLonLatProjection
  module procedure atlas_RotatedLonLatProjection_NP_real64
end interface

!========================================================
contains
!========================================================

function atlas_Projection__ctor_cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Projection) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_Projection__ctor_config(config) result(this)
  use atlas_projection_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Projection) :: this
  type(atlas_Config) :: config
  call this%reset_c_ptr( atlas__Projection__ctor_config(config%CPTR_PGIBUG_B) )
  call this%return()
end function

function atlas_Projection__type(this) result(type_)
  use atlas_projection_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Projection), intent(in) :: this
  character(len=:), allocatable :: type_
  type(c_ptr) :: type_c_str
  integer :: size
  call atlas__Projection__type(this%CPTR_PGIBUG_A, type_c_str, size )
  type_ = c_ptr_to_string(type_c_str)
  call c_ptr_free(type_c_str)
end function

function hash(this)
  use atlas_projection_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Projection), intent(in) :: this
  character(len=:), allocatable :: hash
  type(c_ptr) :: hash_c_str
  integer :: size
  call atlas__Projection__hash(this%CPTR_PGIBUG_A, hash_c_str, size )
  hash = c_ptr_to_string(hash_c_str)
  call c_ptr_free(hash_c_str)
end function

function spec(this)
  use atlas_projection_c_binding
  class(atlas_Projection), intent(in) :: this
  type(atlas_Config) :: spec
  spec = atlas_Config( atlas__Projection__spec(this%CPTR_PGIBUG_A) )
  call spec%return()
end function


!------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Projection__final_auto(this)
  type(atlas_Projection), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Projection__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!------------------------------------------------------------------------------


function atlas_RotatedSchmidtProjection_real64(stretching_factor,north_pole,rotation_angle) result(this)
use, intrinsic :: iso_c_binding, only : c_double
type(atlas_Projection) :: this
real(c_double), intent(in) :: stretching_factor
real(c_double), intent(in), optional :: north_pole(2)
real(c_double), intent(in), optional :: rotation_angle
type(atlas_Config) :: config
config = atlas_Config()
call config%set("type","rotated_schmidt")
call config%set("stretching_factor",stretching_factor)
if( present(rotation_angle) ) then
  call config%set("rotation_angle", rotation_angle)
endif
if( present(north_pole) ) then
  call config%set("north_pole",north_pole)
endif
this = atlas_Projection(config)
call this%return()
end function

! -----------------------------------------------------------------------------

function atlas_LambertConformalConicProjection_real64(longitude0,latitude0,latitude1,latitude2) result(this)
use, intrinsic :: iso_c_binding, only : c_double
type(atlas_Projection) :: this
real(c_double), intent(in) :: longitude0
real(c_double), intent(in) :: latitude0
real(c_double), intent(in), optional :: latitude1
real(c_double), intent(in), optional :: latitude2
type(atlas_Config) :: config
config = atlas_Config()
call config%set("type","lambert_conformal_conic")
call config%set("longitude0",longitude0)
call config%set("latitude0",latitude0)
if( present(latitude1) ) then
  call config%set("latitude1", latitude1)
endif
if( present(latitude2) ) then
  call config%set("latitude2",latitude2)
endif
this = atlas_Projection(config)
call this%return()
end function

! -----------------------------------------------------------------------------

function atlas_RotatedLonLatProjection_NP_real64(north_pole,rotation_angle) result(this)
use, intrinsic :: iso_c_binding, only : c_double
type(atlas_Projection) :: this
real(c_double), intent(in) :: north_pole(2)
real(c_double), intent(in), optional :: rotation_angle
type(atlas_Config) :: config
config = atlas_Config()
call config%set("type","rotated_lonlat")
call config%set("north_pole",north_pole)
if( present(rotation_angle) ) then
  call config%set("rotation_angle", rotation_angle)
endif
this = atlas_Projection(config)
call this%return()
end function

! -----------------------------------------------------------------------------

end module atlas_Projection_module

