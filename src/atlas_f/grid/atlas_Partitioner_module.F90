! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Partitioner_module

use fckit_owned_object_module, only: fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_Partitioner
public :: atlas_MatchingPartitioner
public :: atlas_MatchingMeshPartitioner ! use MatchingPartitioner instead!

private

!-----------------------------
! atlas_Partitioner          !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Partitioner

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
  procedure, public :: partition

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Partitioner__final_auto
#endif

END TYPE atlas_Partitioner

!------------------------------------------------------------------------------

interface atlas_Partitioner
  module procedure atlas_Partitioner__ctor
  module procedure atlas_Partitioner__ctor_type
end interface

interface atlas_MatchingPartitioner
  module procedure atlas_MatchingMeshPartitioner__ctor
  module procedure atlas_MatchingFunctionSpacePartitioner__ctor
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
  call this%reset_c_ptr( atlas__grid__Partitioner__new( config%CPTR_PGIBUG_B ) )
  call this%return()
end function

function atlas_Partitioner__ctor_type( type ) result(this)
  use atlas_config_module, only : atlas_Config
  use atlas_partitioner_c_binding
  use fckit_C_interop_module, only : c_str
  type(atlas_Partitioner) :: this
  character(len=*), intent(in) :: type
  !--------------------------------------
  call this%reset_c_ptr( atlas__grid__Partitioner__new_type( c_str(type) ) )
  call this%return()
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
    call this%reset_c_ptr( atlas__grid__MatchingMeshPartitioner__new( &
      mesh%CPTR_PGIBUG_A, config%CPTR_PGIBUG_B ) )
  else
    opt_config = atlas_Config()
    call this%reset_c_ptr( atlas__grid__MatchingMeshPartitioner__new( &
      mesh%CPTR_PGIBUG_A, opt_config%CPTR_PGIBUG_B ) )
    call opt_config%final()
  endif
  call this%return()
end function


function atlas_MatchingFunctionSpacePartitioner__ctor( functionspace, config ) result(this)
  use atlas_functionspace_module, only : atlas_FunctionSpace
  use atlas_config_module, only : atlas_Config
  use atlas_partitioner_c_binding
  type(atlas_Partitioner) :: this
  class(atlas_FunctionSpace)  , intent(in) :: functionspace
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call this%reset_c_ptr( atlas__grid__MatchingFunctionSpacePartitioner__new( &
      functionspace%CPTR_PGIBUG_A, config%CPTR_PGIBUG_B ) )
  else
    opt_config = atlas_Config()
    call this%reset_c_ptr( atlas__grid__MatchingFunctionSpacePartitioner__new( &
      functionspace%CPTR_PGIBUG_A, opt_config%CPTR_PGIBUG_B ) )
    call opt_config%final()
  endif
  call this%return()
end function

function partition(this,grid) result(distribution)
  use atlas_partitioner_c_binding
  use atlas_GridDistribution_module, only : atlas_GridDistribution
  use atlas_Grid_module, only : atlas_Grid
  type(atlas_GridDistribution) :: distribution
  class(atlas_Partitioner), intent(in) :: this
  class(atlas_Grid), intent(in) :: grid
  distribution = atlas_GridDistribution( atlas__grid__Partitioner__partition( this%CPTR_PGIBUG_A, grid%CPTR_PGIBUG_A ) )
  call distribution%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Partitioner__final_auto(this)
  type(atlas_Partitioner), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Partitioner__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! ----------------------------------------------------------------------------------------

end module atlas_Partitioner_module
