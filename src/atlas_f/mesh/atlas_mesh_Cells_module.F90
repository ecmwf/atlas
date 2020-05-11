! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_mesh_cells_module

use, intrinsic :: iso_c_binding, only: c_ptr
use atlas_HybridElements_module, only: atlas_HybridElements

implicit none

private :: atlas_HybridElements
private :: c_ptr

public :: atlas_mesh_cells

private

type, extends(atlas_HybridElements) :: atlas_mesh_cells
contains
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_mesh_Cells__final_auto
#endif
end type

interface atlas_mesh_cells
  module procedure atlas_mesh_cells__cptr
  module procedure atlas_mesh_cells__constructor
end interface

!========================================================
contains
!========================================================

function atlas_mesh_cells__cptr(cptr) result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_cells) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_mesh_cells__constructor() result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_cells) :: this
  call this%reset_c_ptr( atlas__mesh__HybridElements__create() )
  call this%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_mesh_Cells__final_auto(this)
  type(atlas_mesh_Cells), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_mesh_Cells__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! ----------------------------------------------------------------------------------------

end module atlas_mesh_cells_module

