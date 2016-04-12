
module atlas_mesh_Edges_module

use atlas_HybridElements_module, only: atlas_HybridElements

implicit none

private :: atlas_HybridElements

public :: atlas_mesh_Edges

private

type, extends(atlas_HybridElements) :: atlas_mesh_Edges
contains
end type

interface atlas_mesh_Edges
  module procedure atlas_mesh_Edges__cptr
  module procedure atlas_mesh_Edges__constructor
end interface

!========================================================
contains
!========================================================

function atlas_mesh_edges__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  use atlas_hybridelements_c_binding
  type(atlas_mesh_Edges) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function

function atlas_mesh_edges__constructor() result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_Edges) :: this
  call this%reset_c_ptr( atlas__mesh__HybridElements__create() )
end function

! ----------------------------------------------------------------------------------------

end module atlas_mesh_Edges_module

