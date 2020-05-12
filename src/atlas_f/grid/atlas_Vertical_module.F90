! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Vertical_module

use fckit_shared_object_module, only : fckit_shared_object, fckit_c_deleter, fckit_c_nodeleter
use atlas_field_module, only : atlas_Field

use, intrinsic :: iso_c_binding, only : c_ptr, c_double

implicit none

public :: atlas_Vertical

private

!-----------------------------
! atlas_Vertical     !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(fckit_shared_object) :: atlas_Vertical

! Purpose :
! -------
!   *Vertical* : Object describing vertical levels

! Methods :
! -------

! Author :
! ------
!   28-Nov-2018 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: z
  procedure, public :: size => vsize

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Vertical__final_auto
#endif
END TYPE atlas_Vertical

!------------------------------------------------------------------------------

interface atlas_Vertical
  module procedure atlas_Vertical__ctor_from_cptr
  module procedure atlas_Vertical__ctor_from_array
  module procedure atlas_Vertical__ctor_from_array_with_interval
end interface

private :: c_ptr, c_double
private :: fckit_shared_object
private :: fckit_c_deleter
private :: fckit_c_nodeleter
private :: atlas_Field

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! Vertical routines

function atlas_Vertical__ctor_from_cptr( cptr, delete ) result(this)
  use atlas_vertical_c_binding
  type(c_ptr), value :: cptr
  type(atlas_Vertical) :: this
  logical, optional :: delete
  logical :: opt_delete
  opt_delete = .true.
  if( present(delete) ) opt_delete = delete
  if( opt_delete ) then
    call this%reset_c_ptr( cptr, fckit_c_deleter(atlas__Vertical__delete) )
  else
    call this%reset_c_ptr( cptr, fckit_c_nodeleter() )
  endif
  call this%return()
end function

function atlas_Vertical__ctor_from_array( levels ) result(this)
  use atlas_vertical_c_binding
  use atlas_kinds_module, only : ATLAS_KIND_IDX
  type(atlas_Vertical) :: this
  real(c_double), intent(in) :: levels(:)
  integer(ATLAS_KIND_IDX) :: nb_levels
  nb_levels = size(levels)
  call this%reset_c_ptr( atlas__Vertical__new(nb_levels, levels), &
    & fckit_c_deleter(atlas__Vertical__delete) )
  call this%return()
end function

function atlas_Vertical__ctor_from_array_with_interval( levels, interval ) result(this)
  use atlas_vertical_c_binding
  use atlas_kinds_module, only : ATLAS_KIND_IDX
  type(atlas_Vertical) :: this
  real(c_double), intent(in) :: levels(:)
  real(c_double), intent(in) :: interval(2)
  integer(ATLAS_KIND_IDX) :: nb_levels
  nb_levels = size(levels)
  call this%reset_c_ptr( atlas__Vertical__new_interval(nb_levels, levels, interval), &
    & fckit_c_deleter(atlas__Vertical__delete) )
  call this%return()
end function

! ----------------------------------------------------------------------------------------

function z( this )
  use atlas_vertical_c_binding
  type(atlas_Field) :: z
  class(atlas_Vertical), intent(in) :: this
  z = atlas_Field( atlas__Vertical__z( this%CPTR_PGIBUG_B ) )
  call z%return()
end function

! ----------------------------------------------------------------------------------------

function vsize( this )
  use atlas_vertical_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  integer(c_int) :: vsize
  class(atlas_Vertical), intent(in) :: this
  vsize = atlas__Vertical__size( this%CPTR_PGIBUG_B )
end function

! ----------------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Vertical__final_auto(this)
  type(atlas_Vertical), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Vertical__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! ----------------------------------------------------------------------------------------

end module atlas_Vertical_module
