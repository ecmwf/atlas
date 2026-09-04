! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Mesh_module

use fckit_owned_object_module, only: fckit_owned_object
use atlas_mesh_Cells_module, only: atlas_mesh_Cells
use atlas_mesh_Nodes_module, only: atlas_mesh_Nodes
use atlas_mesh_Edges_module, only: atlas_mesh_Edges
use atlas_Grid_module, only: atlas_Grid
use atlas_GridDistribution_module, only: atlas_GridDistribution
use atlas_Partitioner_module, only: atlas_Partitioner
use atlas_Config_module, only: atlas_Config
use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr

implicit none

private :: fckit_owned_object
private :: atlas_mesh_Cells
private :: atlas_mesh_Nodes
private :: atlas_mesh_Edges
private :: atlas_GridDistribution
private :: atlas_Partitioner
private :: atlas_Config
private :: c_size_t
private :: c_ptr

public :: atlas_Mesh

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

TYPE, extends(fckit_owned_object) :: atlas_Mesh

! Purpose :
! -------
!   *Mesh* : Container type holding an entire mesh

! Methods :
! -------

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, public :: nodes => Mesh__nodes
  procedure, public :: cells => Mesh__cells
  procedure, public :: edges => Mesh__edges
  procedure, public :: footprint

  procedure, public :: update_device
  procedure, public :: update_host
  procedure, public :: sync_host_device

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Mesh__final_auto
#endif
END TYPE atlas_Mesh

interface atlas_Mesh
  module procedure atlas_Mesh__cptr
  module procedure atlas_Mesh__ctor
  module procedure atlas_Mesh__ctor_grid
  module procedure atlas_Mesh__ctor_grid_distribution
  module procedure atlas_Mesh__ctor_grid_partitioner
end interface

!========================================================
contains
!========================================================

function atlas_Mesh__cptr(cptr) result(this)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function atlas_Mesh__cptr

!-------------------------------------------------------------------------------

function atlas_Mesh__ctor() result(this)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: this
  call this%reset_c_ptr( atlas__Mesh__new() )
  call this%return()
end function atlas_Mesh__ctor

!-------------------------------------------------------------------------------

function atlas_Mesh__ctor_grid(grid) result(this)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: this
  class(atlas_Grid), intent(in) :: grid
  call this%reset_c_ptr( atlas__Mesh__new_grid(grid%c_ptr()) )
  call this%return()
end function atlas_Mesh__ctor_grid

!-------------------------------------------------------------------------------

function atlas_Mesh__ctor_grid_distribution(grid, distribution, config) result(this)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: this
  class(atlas_Grid), intent(in) :: grid
  class(atlas_GridDistribution), intent(in) :: distribution
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if (present(config)) then
    call this%reset_c_ptr( atlas__Mesh__new_grid_distribution( &
      grid%c_ptr(), distribution%c_ptr(), config%c_ptr()) )
  else
    opt_config = atlas_Config()
    call this%reset_c_ptr( atlas__Mesh__new_grid_distribution( &
      grid%c_ptr(), distribution%c_ptr(), opt_config%c_ptr()) )
    call opt_config%final()
  endif
  call this%return()
end function atlas_Mesh__ctor_grid_distribution

!-------------------------------------------------------------------------------

function atlas_Mesh__ctor_grid_partitioner(grid, partitioner, config) result(this)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: this
  class(atlas_Grid), intent(in) :: grid
  class(atlas_Partitioner), intent(in) :: partitioner
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if (present(config)) then
    call this%reset_c_ptr( atlas__Mesh__new_grid_partitioner( &
      grid%c_ptr(), partitioner%c_ptr(), config%c_ptr()) )
  else
    opt_config = atlas_Config()
    call this%reset_c_ptr( atlas__Mesh__new_grid_partitioner( &
      grid%c_ptr(), partitioner%c_ptr(), opt_config%c_ptr()) )
    call opt_config%final()
  endif
  call this%return()
end function atlas_Mesh__ctor_grid_partitioner

!-------------------------------------------------------------------------------

function Mesh__nodes(this) result(nodes)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Nodes) :: nodes
  nodes = atlas_mesh_Nodes( atlas__Mesh__nodes(this%c_ptr()) )
  call nodes%return()
end function

!-------------------------------------------------------------------------------

function Mesh__cells(this) result(cells)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Cells) :: cells
  cells = atlas_mesh_Cells( atlas__Mesh__cells(this%c_ptr()) )
  call cells%return()
end function

!-------------------------------------------------------------------------------

function Mesh__edges(this) result(edges)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Edges) :: edges
  edges = atlas_mesh_Edges( atlas__Mesh__Edges(this%c_ptr()) )
  call edges%return()
end function

!-------------------------------------------------------------------------------

function footprint(this)
  use atlas_mesh_c_binding
  integer(c_size_t) :: footprint
  class(atlas_Mesh) :: this
  footprint = atlas__Mesh__footprint(this%c_ptr())
end function

!-------------------------------------------------------------------------------

subroutine update_device(this)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(inout) :: this
  call atlas__Mesh__update_device(this%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

subroutine update_host(this)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(inout) :: this
  call atlas__Mesh__update_host(this%c_ptr())
end subroutine

! ----------------------------------------------------------------------------------------

subroutine sync_host_device(this)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(inout) :: this
  call atlas__Mesh__sync_host_device(this%c_ptr())
end subroutine

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Mesh__final_auto(this)
  type(atlas_Mesh), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Mesh__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!-------------------------------------------------------------------------------

end module atlas_Mesh_module
