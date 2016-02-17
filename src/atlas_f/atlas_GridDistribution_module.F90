
module atlas_GridDistribution_module


use atlas_object_module, only: atlas_object

implicit none

private :: atlas_object

public :: atlas_GridDistribution

private

!-----------------------------
! atlas_GridDistribution     !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_GridDistribution

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

END TYPE atlas_GridDistribution

!------------------------------------------------------------------------------

interface atlas_GridDistribution
  module procedure atlas_GridDistribution__ctor
end interface

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! GridDistribution routines

function atlas_GridDistribution__ctor( part, part0 ) result(griddistribution)
  use atlas_griddistribution_c_binding
  type(atlas_GridDistribution) :: griddistribution
  integer, intent(in) :: part(:)
  integer, intent(in), optional :: part0
  integer:: npts, opt_part0
  opt_part0 = 0
  if( present(part0) ) opt_part0 = part0
  npts = size(part)
  call griddistribution%reset_c_ptr( atlas__GridDistribution__new(npts, part, opt_part0) );
end function atlas_GridDistribution__ctor


subroutine atlas_GridDistribution__final( this )
  use atlas_griddistribution_c_binding
  class(atlas_GridDistribution), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__GridDistribution__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_GridDistribution__delete( this )
  use atlas_griddistribution_c_binding
  class(atlas_GridDistribution), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__GridDistribution__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_GridDistribution_module
