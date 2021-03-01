! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Domain_module

use, intrinsic :: iso_c_binding, only : c_double
use fckit_owned_object_module, only : fckit_owned_object
use atlas_config_module, only : atlas_Config

implicit none

private :: fckit_owned_object
private :: atlas_Config

public :: atlas_Domain
public :: atlas_LonLatRectangularDomain

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Domain

! Purpose :
! -------
!   *atlas_Domain* : 

!------------------------------------------------------------------------------
contains
  procedure, public :: type => atlas_Domain__type
  procedure, public :: hash
  procedure, public :: spec

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Domain__final_auto
#endif

END TYPE atlas_Domain

interface atlas_Domain
  module procedure atlas_Domain__ctor_cptr
  module procedure atlas_Domain__ctor_config
end interface

!------------------------------------------------------------------------------
TYPE, extends(atlas_Domain) :: atlas_LonLatRectangularDomain

! Purpose :
! -------
!   *atlas_LonLatRectangularDomain* : 

!------------------------------------------------------------------------------
contains
  procedure, public :: north
  procedure, public :: west
  procedure, public :: south
  procedure, public :: east
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_LonLatRectangularDomain__final_auto
#endif

END TYPE atlas_LonLatRectangularDomain

interface atlas_LonLatRectangularDomain
  module procedure atlas_LonLatRectangularDomain__ctor_cptr
end interface
!------------------------------------------------------------------------------


!========================================================
contains
!========================================================

function atlas_Domain__ctor_cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Domain) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_Domain__ctor_config(config) result(this)
  use atlas_Domain_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Domain) :: this
  type(atlas_Config) :: config
  call this%reset_c_ptr( atlas__Domain__ctor_config(config%CPTR_PGIBUG_B) )
  call this%return()
end function

function atlas_Domain__type(this) result(type_)
  use atlas_Domain_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Domain), intent(in) :: this
  character(len=:), allocatable :: type_
  type(c_ptr) :: type_c_str
  integer :: size
  call atlas__Domain__type(this%CPTR_PGIBUG_A, type_c_str, size )
  type_ = c_ptr_to_string(type_c_str)
  call c_ptr_free(type_c_str)
end function

function hash(this)
  use atlas_Domain_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Domain), intent(in) :: this
  character(len=:), allocatable :: hash
  type(c_ptr) :: hash_c_str
  integer :: size
  call atlas__Domain__hash(this%CPTR_PGIBUG_A, hash_c_str, size )
  hash = c_ptr_to_string(hash_c_str)
  call c_ptr_free(hash_c_str)
end function

function spec(this)
  use atlas_Domain_c_binding
  class(atlas_Domain), intent(in) :: this
  type(atlas_Config) :: spec
  spec = atlas_Config( atlas__Domain__spec(this%CPTR_PGIBUG_A) )
  call spec%return()
end function


!------------------------------------------------------------------------------

function atlas_LonLatRectangularDomain__ctor_cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_LonLatRectangularDomain) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function north(this) result(value)
  use atlas_Domain_c_binding
  class(atlas_LonLatRectangularDomain), intent(in) :: this
  real(c_double) :: value
  value = atlas__LonLatRectangularDomain__north(this%CPTR_PGIBUG_A)
end function

function west(this) result(value)
  use atlas_Domain_c_binding
  class(atlas_LonLatRectangularDomain), intent(in) :: this
  real(c_double) :: value
  value = atlas__LonLatRectangularDomain__west(this%CPTR_PGIBUG_A)
end function

function south(this) result(value)
  use atlas_Domain_c_binding
  class(atlas_LonLatRectangularDomain), intent(in) :: this
  real(c_double) :: value
  value = atlas__LonLatRectangularDomain__south(this%CPTR_PGIBUG_A)
end function

function east(this) result(value)
  use atlas_Domain_c_binding
  class(atlas_LonLatRectangularDomain), intent(in) :: this
  real(c_double) :: value
  value = atlas__LonLatRectangularDomain__east(this%CPTR_PGIBUG_A)
end function

!------------------------------------------------------------------------------


#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Domain__final_auto(this)
  type(atlas_Domain), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Domain__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_LonLatRectangularDomain__final_auto(this)
  type(atlas_LonLatRectangularDomain), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_LonLatRectangularDomain__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif


end module atlas_Domain_module

