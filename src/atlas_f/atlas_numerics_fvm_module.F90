
module atlas_numerics_fvm_module

use atlas_refcounted_module, only : atlas_RefCounted
use atlas_numerics_Method_module, only : atlas_numerics_Method
use atlas_Config_module, only : atlas_Config
use atlas_Mesh_module, only : atlas_Mesh
use atlas_functionspace_Nodes_module, only : atlas_functionspace_Nodes
use atlas_functionspace_Edges_module, only : atlas_functionspace_Edges
implicit none

private :: atlas_RefCounted
private :: atlas_numerics_Method
private :: atlas_Mesh
private :: atlas_Config
private :: atlas_functionspace_Nodes

public :: atlas_numerics_fvm_Method

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_numerics_Method) :: atlas_numerics_fvm_Method

! Purpose :
! -------
!   *Method* :
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done

! Methods :
! -------
!   name : The name or tag this function space was created with

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, public :: functionspace_nodes
  procedure, public :: functionspace_edges

END TYPE atlas_numerics_fvm_Method

interface atlas_numerics_fvm_Method
  module procedure atlas_numerics_fvm_Method__cptr
  module procedure atlas_numerics_fvm_Method__mesh_config
end interface

!========================================================
contains
!========================================================

function atlas_numerics_fvm_Method__cptr(cptr) result(method)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_numerics_fvm_Method) :: method
  type(c_ptr), intent(in) :: cptr
  call method%reset_c_ptr( cptr )
end function

function atlas_numerics_fvm_Method__mesh_config(mesh,config) result(method)
  use atlas_numerics_fvm_method_c_binding
  type(atlas_numerics_fvm_Method) :: method
  type(atlas_Mesh), intent(inout) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    method = atlas_numerics_fvm_Method__cptr( &
      & atlas__numerics__fvm__Method__new(mesh%c_ptr(),config%c_ptr()) )
  else
    opt_config = atlas_Config()
    method = atlas_numerics_fvm_Method__cptr( &
      & atlas__numerics__fvm__Method__new(mesh%c_ptr(),opt_config%c_ptr()) )
    call opt_config%final()
  endif
  call method%return()
end function

function functionspace_nodes(this)
  use atlas_numerics_fvm_method_c_binding
  type(atlas_functionspace_Nodes) :: functionspace_nodes
  class(atlas_numerics_fvm_Method) :: this
  functionspace_nodes = atlas_functionspace_Nodes( &
    & atlas__numerics__fvm__Method__functionspace_nodes(this%c_ptr()) )
  call functionspace_nodes%return()
end function

function functionspace_edges(this)
  use atlas_numerics_fvm_method_c_binding
  type(atlas_functionspace_Nodes) :: functionspace_edges
  class(atlas_numerics_fvm_Method) :: this
  functionspace_edges = atlas_functionspace_Edges( &
    & atlas__numerics__fvm__Method__functionspace_edges(this%c_ptr()) )
  call functionspace_edges%return()
end function


end module atlas_numerics_fvm_module

