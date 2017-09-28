
module atlas_Partitioner_module

use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_Partitioner
public :: atlas_MatchingMeshPartitioner

private

!-----------------------------
! atlas_Partitioner          !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_Partitioner

! Purpose :
! -------
!   *Partitioner* : Object passed from atlas to inspect grid distribution

! Methods :
! -------

! Author :
! ------
!   12-Mar-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, public :: final  => atlas_Partitioner__final
  procedure, public :: delete => atlas_Partitioner__delete
  procedure, public :: copy   => atlas_Partitioner__copy

  procedure, public :: partition

END TYPE atlas_Partitioner

!------------------------------------------------------------------------------

interface atlas_Partitioner
  module procedure atlas_Partitioner__ctor
end interface

interface atlas_MatchingMeshPartitioner
  module procedure atlas_MatchingMeshPartitioner__ctor
end interface

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! Partitioner routines

function atlas_Partitioner__ctor( config ) result(this)
  use atlas_config_module, only : atlas_Config
  use atlas_partitioner_c_binding
  type(atlas_Partitioner) :: this
  type(atlas_Config) :: config
  call this%reset_c_ptr( atlas__grid__Partitioner__new( config%c_ptr() ) )
end function

function atlas_MatchingMeshPartitioner__ctor( mesh, config ) result(this)
  use atlas_mesh_module, only : atlas_Mesh
  use atlas_config_module, only : atlas_Config
  use atlas_partitioner_c_binding
  type(atlas_Partitioner) :: this
  type(atlas_Mesh)  , intent(in) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call this%reset_c_ptr( atlas__grid__MatchingMeshPartitioner__new( mesh%c_ptr(), config%c_ptr() ) )
   else
     opt_config = atlas_Config()
     call this%reset_c_ptr( atlas__grid__MatchingMeshPartitioner__new( mesh%c_ptr(), opt_config%c_ptr() ) )
     call opt_config%final()
   endif
end function

subroutine atlas_Partitioner__final( this )
  use atlas_partitioner_c_binding
  class(atlas_Partitioner), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__grid__Partitioner__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_Partitioner__delete( this )
  use atlas_partitioner_c_binding
  class(atlas_Partitioner), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__grid__Partitioner__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_Partitioner__copy(this,obj_in)
  class(atlas_Partitioner), intent(inout) :: this
  class(fckit_RefCounted), target, intent(in) :: obj_in
end subroutine

function partition(this,grid) result(distribution)
  use atlas_partitioner_c_binding
  use atlas_GridDistribution_module, only : atlas_GridDistribution
  use atlas_Grid_module, only : atlas_Grid
  type(atlas_GridDistribution) :: distribution
  class(atlas_Partitioner), intent(in) :: this
  class(atlas_Grid), intent(in) :: grid
  distribution = atlas_GridDistribution( atlas__grid__Partitioner__partition( this%c_ptr(), grid%c_ptr() ) )
  call distribution%return()
end function

! ----------------------------------------------------------------------------------------

end module atlas_Partitioner_module
