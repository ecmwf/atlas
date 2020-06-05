! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_fvm_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_Method_module, only : atlas_Method
implicit none

private :: fckit_owned_object
private :: atlas_Method

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

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_fvm_Method__final_auto
#endif

END TYPE atlas_fvm_Method

interface atlas_fvm_Method
  module procedure atlas_fvm_Method__cptr
  module procedure atlas_fvm_Method__mesh_config
end interface

!========================================================
contains
!========================================================

function atlas_fvm_Method__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_fvm_Method) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function

function atlas_fvm_Method__mesh_config(mesh,config) result(this)
  use atlas_fvm_method_c_binding
  use atlas_Config_module, only : atlas_Config
  use atlas_Mesh_module, only : atlas_Mesh
  type(atlas_fvm_Method) :: this
  type(atlas_Mesh), intent(inout) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call this%reset_c_ptr( atlas__numerics__fvm__Method__new(mesh%CPTR_PGIBUG_A, &
      config%CPTR_PGIBUG_B) )
  else
    opt_config = atlas_Config()
    call this%reset_c_ptr( atlas__numerics__fvm__Method__new(mesh%CPTR_PGIBUG_A, &
      opt_config%CPTR_PGIBUG_B) )
    call opt_config%final()
  endif
  call this%return()
end function

function node_columns(this)
  use atlas_fvm_method_c_binding
  use atlas_functionspace_NodeColumns_module, only : atlas_functionspace_NodeColumns
  type(atlas_functionspace_NodeColumns) :: node_columns
  class(atlas_fvm_Method) :: this
  node_columns = atlas_functionspace_NodeColumns( &
    & atlas__numerics__fvm__Method__functionspace_nodes(this%CPTR_PGIBUG_A) )
  call node_columns%return()
end function

function edge_columns(this)
  use atlas_fvm_method_c_binding
  use atlas_functionspace_EdgeColumns_module, only : atlas_functionspace_EdgeColumns
  type(atlas_functionspace_EdgeColumns) :: edge_columns
  class(atlas_fvm_Method) :: this
  edge_columns = atlas_functionspace_EdgeColumns( &
    & atlas__numerics__fvm__Method__functionspace_edges(this%CPTR_PGIBUG_A) )
  call edge_columns%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_fvm_Method__final_auto(this)
  type(atlas_fvm_Method), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_fvm_Method__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_fvm_module

