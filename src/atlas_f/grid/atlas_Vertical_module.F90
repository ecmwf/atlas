#include "atlas/atlas_f.h"

module atlas_Vertical_module

use fckit_shared_object_module, only : fckit_shared_object, fckit_c_deleter, fckit_c_nodeleter

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
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Vertical__final_auto
#endif
END TYPE atlas_Vertical

!------------------------------------------------------------------------------

interface atlas_Vertical
  module procedure atlas_Vertical__ctor_from_cptr
  module procedure atlas_Vertical__ctor_from_array
end interface

private :: c_ptr, c_double
private :: fckit_shared_object
private :: fckit_c_deleter
private :: fckit_c_nodeleter

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
  real(c_double) :: levels(:)
  integer(ATLAS_KIND_IDX) :: nb_levels
  nb_levels = size(levels)
  call this%reset_c_ptr( atlas__Vertical__new(nb_levels, levels), &
    & fckit_c_deleter(atlas__Vertical__delete) )
  call this%return()
end function

! ----------------------------------------------------------------------------------------

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

! ----------------------------------------------------------------------------------------

end module atlas_Vertical_module
