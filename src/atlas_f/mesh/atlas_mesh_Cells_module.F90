
module atlas_mesh_cells_module

use atlas_HybridElements_module, only: atlas_HybridElements

implicit none

private :: atlas_HybridElements

public :: atlas_mesh_cells

private

type, extends(atlas_HybridElements) :: atlas_mesh_cells
contains
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
  use, intrinsic :: iso_c_binding, only: c_ptr
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

! ----------------------------------------------------------------------------------------

end module atlas_mesh_cells_module

