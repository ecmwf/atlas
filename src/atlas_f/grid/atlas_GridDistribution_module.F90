
module atlas_GridDistribution_module

use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_GridDistribution

private

!-----------------------------
! atlas_GridDistribution     !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_GridDistribution

! Purpose :
! -------
!   *GridDistribution* : Object passed from atlas to inspect grid distribution

! Methods :
! -------

! Author :
! ------
!   12-Mar-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, public :: final  => atlas_GridDistribution__final
  procedure, public :: delete => atlas_GridDistribution__delete
  procedure, public :: copy   => atlas_GridDistribution__copy

END TYPE atlas_GridDistribution

!------------------------------------------------------------------------------

interface atlas_GridDistribution
  module procedure atlas_GridDistribution__cptr
  module procedure atlas_GridDistribution__ctor
end interface

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! GridDistribution routines

function atlas_GridDistribution__cptr( cptr ) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_distribution_c_binding
  type(atlas_GridDistribution) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr );
end function

function atlas_GridDistribution__ctor( part, part0 ) result(this)
  use atlas_distribution_c_binding
  type(atlas_GridDistribution) :: this
  integer, intent(in) :: part(:)
  integer, intent(in), optional :: part0
  integer:: npts, opt_part0
  opt_part0 = 0
  if( present(part0) ) opt_part0 = part0
  npts = size(part)
  this = atlas_GridDistribution__cptr( atlas__GridDistribution__new(npts, part, opt_part0) )
  call this%return()
end function


subroutine atlas_GridDistribution__final( this )
  use atlas_distribution_c_binding
  class(atlas_GridDistribution), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__GridDistribution__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_GridDistribution__delete( this )
  use atlas_distribution_c_binding
  class(atlas_GridDistribution), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__GridDistribution__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_GridDistribution__copy(this,obj_in)
  class(atlas_GridDistribution), intent(inout) :: this
  class(fckit_RefCounted), target, intent(in) :: obj_in
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_GridDistribution_module
