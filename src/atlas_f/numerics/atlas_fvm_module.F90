
module atlas_fvm_module

use fckit_refcounted_module, only : fckit_refcounted
use atlas_Method_module, only : atlas_Method
use atlas_Config_module, only : atlas_Config
use atlas_Mesh_module, only : atlas_Mesh
use atlas_functionspace_NodeColumns_module, only : atlas_functionspace_NodeColumns
use atlas_functionspace_EdgeColumns_module, only : atlas_functionspace_EdgeColumns
implicit none

private :: fckit_refcounted
private :: atlas_Method
private :: atlas_Mesh
private :: atlas_Config
private :: atlas_functionspace_NodeColumns

public :: atlas_fvm_Method

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_Method) :: atlas_fvm_Method

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

  procedure, public :: node_columns
  procedure, public :: edge_columns

END TYPE atlas_fvm_Method

interface atlas_fvm_Method
  module procedure atlas_fvm_Method__cptr
  module procedure atlas_fvm_Method__mesh_config
end interface

!========================================================
contains
!========================================================

function atlas_fvm_Method__cptr(cptr) result(method)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_fvm_Method) :: method
  type(c_ptr), intent(in) :: cptr
  call method%reset_c_ptr( cptr )
end function

function atlas_fvm_Method__mesh_config(mesh,config) result(method)
  use atlas_fvm_method_c_binding
  type(atlas_fvm_Method) :: method
  type(atlas_Mesh), intent(inout) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    method = atlas_fvm_Method__cptr( &
      & atlas__numerics__fvm__Method__new(mesh%c_ptr(),config%c_ptr()) )
  else
    opt_config = atlas_Config()
    method = atlas_fvm_Method__cptr( &
      & atlas__numerics__fvm__Method__new(mesh%c_ptr(),opt_config%c_ptr()) )
    call opt_config%final()
  endif
  call method%return()
end function

function node_columns(this)
  use atlas_fvm_method_c_binding
  type(atlas_functionspace_NodeColumns) :: node_columns
  class(atlas_fvm_Method) :: this
  node_columns = atlas_functionspace_NodeColumns( &
    & atlas__numerics__fvm__Method__functionspace_nodes(this%c_ptr()) )
  call node_columns%return()
end function

function edge_columns(this)
  use atlas_fvm_method_c_binding
  type(atlas_functionspace_NodeColumns) :: edge_columns
  class(atlas_fvm_Method) :: this
  edge_columns = atlas_functionspace_EdgeColumns( &
    & atlas__numerics__fvm__Method__functionspace_edges(this%c_ptr()) )
  call edge_columns%return()
end function


end module atlas_fvm_module

