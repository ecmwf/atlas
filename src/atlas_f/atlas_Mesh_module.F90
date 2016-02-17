
#include "atlas_f/atlas_f_defines.h"

module atlas_Mesh_module


use, intrinsic :: iso_c_binding, only: c_ptr
use atlas_c_interop, only: c_str
use atlas_refcounted_module, only: atlas_refcounted
use atlas_mesh_HybridElements_module, only: atlas_mesh_Cells, atlas_mesh_Edges
use atlas_mesh_Nodes_module, only: atlas_mesh_Nodes

#if ! DEPRECATE_OLD_FUNCTIONSPACE
use atlas_deprecated_FunctionSpace_module, only : atlas_deprecated_FunctionSpace, &
    & ATLAS_FIELD_NB_VARS
#endif

implicit none

private :: c_ptr
private :: c_str
private :: atlas_refcounted
private :: atlas_mesh_Cells, atlas_mesh_Edges
private :: atlas_mesh_Nodes

public :: atlas_Mesh

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

TYPE, extends(atlas_RefCounted) :: atlas_Mesh

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
  procedure :: create_nodes => Mesh__create_nodes
  procedure :: nodes => Mesh__nodes
  procedure :: cells => Mesh__cells
  procedure :: edges => Mesh__edges
  procedure, public :: delete => atlas_Mesh__delete
  procedure, public :: copy => atlas_Mesh__copy

#if ! DEPRECATE_OLD_FUNCTIONSPACE
  procedure, private :: create_function_space_nodes => Mesh__create_function_space_nodes
  procedure, private :: create_function_space_shape => Mesh__create_function_space_shape
  generic :: create_function_space => create_function_space_nodes, create_function_space_shape
  procedure :: function_space => Mesh__function_space
#endif

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

#if !DEPRECATE_OLD_FUNCTIONSPACE
subroutine Mesh__create_function_space_nodes(this,name,shape_func,nb_nodes)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_nodes
  integer :: shape(2)
  integer , parameter :: fortran_ordering = 1
  shape = (/ATLAS_FIELD_NB_VARS,nb_nodes/)
  call atlas__Mesh__create_function_space(this%c_ptr(),c_str(name),c_str(shape_func), &
  & shape, size(shape), fortran_ordering )
end subroutine Mesh__create_function_space_nodes

subroutine Mesh__create_function_space_shape(this,name,shape_func,shape)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: shape(:)
  integer , parameter :: fortran_ordering = 1
  call atlas__Mesh__create_function_space(this%c_ptr(),c_str(name),c_str(shape_func), &
  & shape, size(shape), fortran_ordering )
end subroutine Mesh__create_function_space_shape

function Mesh__function_space(this,name) result(function_space)
  use atlas_mesh_c_binding
  class(atlas_Mesh), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_deprecated_FunctionSpace) :: function_space
  call function_space%reset_c_ptr( atlas__Mesh__function_space(this%c_ptr(), c_str(name) ) )
  if( function_space%is_null() ) write(0,*) 'call abort()'
end function Mesh__function_space
#endif

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
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_Mesh_module
