
#include "atlas/atlas_f.h"

module atlas_Mesh_module


use, intrinsic :: iso_c_binding, only: c_ptr
use fckit_c_interop_module, only: c_str
use fckit_refcounted_module, only: fckit_refcounted
use atlas_mesh_Cells_module, only: atlas_mesh_Cells
use atlas_mesh_Edges_module, only: atlas_mesh_Edges
use atlas_mesh_Nodes_module, only: atlas_mesh_Nodes

implicit none

private :: c_ptr
private :: c_str
private :: fckit_refcounted
private :: atlas_mesh_Cells
private :: atlas_mesh_Edges
private :: atlas_mesh_Nodes

public :: atlas_Mesh

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

TYPE, extends(fckit_refcounted) :: atlas_Mesh

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
  procedure, public :: create_nodes => Mesh__create_nodes
  procedure, public :: nodes => Mesh__nodes
  procedure, public :: cells => Mesh__cells
  procedure, public :: edges => Mesh__edges
  procedure, public :: delete => atlas_Mesh__delete
  procedure, public :: copy => atlas_Mesh__copy
  procedure, public :: footprint

END TYPE atlas_Mesh

interface atlas_Mesh
  module procedure atlas_Mesh__cptr
  module procedure atlas_Mesh__ctor
end interface

!========================================================
contains
!========================================================

function atlas_Mesh__cptr(cptr) result(mesh)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: mesh
  type(c_ptr), intent(in) :: cptr
  call mesh%reset_c_ptr( cptr )
end function atlas_Mesh__cptr

function atlas_Mesh__ctor() result(mesh)
  use atlas_mesh_c_binding
  type(atlas_Mesh) :: mesh
  call mesh%reset_c_ptr( atlas__Mesh__new() )
end function atlas_Mesh__ctor

function Mesh__create_nodes(this,nb_nodes) result(nodes)
  use atlas_mesh_c_binding
  type(atlas_mesh_Nodes) :: nodes
  class(atlas_Mesh), intent(in) :: this
  integer, intent(in) :: nb_nodes
  call nodes%reset_c_ptr( atlas__Mesh__create_nodes(this%c_ptr(),nb_nodes) )
  if( nodes%is_null() ) write(0,*) 'call abort()'
end function

function Mesh__nodes(this) result(nodes)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Nodes) :: nodes
  call nodes%reset_c_ptr( atlas__Mesh__nodes(this%c_ptr()) )
  if( nodes%is_null() ) write(0,*) 'call abort()'
end function

function Mesh__cells(this) result(cells)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Cells) :: cells
  cells = atlas_mesh_Cells(atlas__Mesh__cells(this%c_ptr()))
  if( cells%is_null() ) write(0,*) 'call abort()'
  call cells%return()
end function

function Mesh__edges(this) result(edges)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  type(atlas_mesh_Edges) :: edges
  edges = atlas_mesh_Edges( atlas__Mesh__Edges(this%c_ptr()) )
  if( edges%is_null() ) write(0,*) 'call abort()'
  call edges%return()
end function

subroutine atlas_Mesh__delete(this)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Mesh__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_Mesh__delete

subroutine atlas_Mesh__copy(this,obj_in)
  class(atlas_Mesh), intent(inout) :: this
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine

function footprint(this)
  use, intrinsic :: iso_c_binding, only : c_size_t
  use atlas_mesh_c_binding
  integer(c_size_t) :: footprint
  class(atlas_Mesh) :: this
  footprint = atlas__Mesh__footprint(this%c_ptr())
end function

! ----------------------------------------------------------------------------------------

end module atlas_Mesh_module
