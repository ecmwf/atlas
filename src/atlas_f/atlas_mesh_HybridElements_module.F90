
module atlas_mesh_hybridelements_module

!use iso_c_binding, only : c_funptr, c_ptr, c_loc, c_f_pointer, c_f_procpointer, c_funloc, c_int, c_size_t
use atlas_refcounted_module
use atlas_Connectivity_module, only: atlas_MultiBlockConnectivity
use iso_c_binding, only: c_size_t, c_int

implicit none

private :: c_size_t, c_int
private :: atlas_MultiBlockConnectivity

!-----------------------------
! atlas_mesh_HybridElements  !
!-----------------------------

type, extends(atlas_refcounted), public :: atlas_mesh_HybridElements
contains

! Public methods
  procedure, public :: copy     => atlas_mesh_HybridElements__copy
  procedure, public :: delete   => atlas_mesh_HybridElements__delete
  procedure, public :: size     => atlas_mesh_HybridElements__size

  procedure, public ::  node_connectivity => atlas_mesh_HybridElements__node_connectivity
  procedure, public ::  edge_connectivity => atlas_mesh_HybridElements__edge_connectivity
  procedure, public ::  cell_connectivity => atlas_mesh_HybridElements__cell_connectivity

  generic, public :: add => add_elements, add_elements_with_nodes

!  generic, public :: field => field_by_idx, field_by_name

!  procedure, public :: nb_fields
!  procedure, public :: has_field
!  procedure, public :: global_index
!  procedure, public :: remote_index
!  procedure, public :: partition
!  procedure, public :: halo

! Private methods
  procedure, private :: add_elements
  procedure, private :: add_elements_with_nodes
!  procedure, private :: add_field
!  procedure, private :: field_by_idx
!  procedure, private :: field_by_name

end type

interface atlas_mesh_HybridElements
  module procedure atlas_mesh_HybridElements__constructor
end interface

type, extends(atlas_mesh_HybridElements), public :: atlas_mesh_Edges
contains
end type

interface atlas_mesh_Edges
  module procedure atlas_mesh_HybridElements__constructor
end interface

type, extends(atlas_mesh_HybridElements), public :: atlas_mesh_Cells
contains
end type

interface atlas_mesh_Cells
  module procedure atlas_mesh_HybridElements__constructor
end interface

!-----------------------------
! atlas_mesh_ElementType     !
!-----------------------------

type, extends(atlas_refcounted), public :: atlas_mesh_ElementType
contains
! Public methods
  procedure, public :: copy     => atlas_mesh_ElementType__copy
  procedure, public :: delete   => atlas_mesh_ElementType__delete
end type

interface atlas_mesh_Triangle
  module procedure atlas_mesh_Triangle__constructor
end interface

interface atlas_mesh_Quadrilateral
  module procedure atlas_mesh_Quadrilateral__constructor
end interface

interface atlas_mesh_Line
  module procedure atlas_mesh_Line__constructor
end interface


!========================================================
contains
!========================================================

function atlas_mesh_HybridElements__constructor() result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_HybridElements) :: this
  call this%reset_c_ptr( atlas__mesh__HybridElements__create() )
end function

subroutine atlas_mesh_HybridElements__delete(this)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_HybridElements), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__mesh__HybridElements__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_mesh_HybridElements__copy(this,obj_in)
  class(atlas_mesh_HybridElements), intent(inout) :: this
  class(atlas_RefCounted),   target, intent(in) :: obj_in
end subroutine

function atlas_mesh_HybridElements__size(this) result(val)
  use atlas_hybridelements_c_binding
  integer(c_size_t) :: val
  class(atlas_mesh_HybridElements), intent(in) :: this
  val = atlas__mesh__HybridElements__size(this%c_ptr())
end function

function atlas_mesh_HybridElements__node_connectivity(this) result(connectivity)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_HybridElements), intent(in) :: this
  type(atlas_MultiBlockConnectivity) :: connectivity
  connectivity = atlas_MultiBlockConnectivity( &
      atlas__mesh__HybridElements__node_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function atlas_mesh_HybridElements__edge_connectivity(this) result(connectivity)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_HybridElements), intent(in) :: this
  type(atlas_MultiBlockConnectivity) :: connectivity
  connectivity = atlas_MultiBlockConnectivity( &
      atlas__mesh__HybridElements__edge_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function atlas_mesh_HybridElements__cell_connectivity(this) result(connectivity)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_HybridElements), intent(in) :: this
  type(atlas_MultiBlockConnectivity) :: connectivity
  connectivity = atlas_MultiBlockConnectivity( &
      atlas__mesh__HybridElements__cell_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

subroutine add_elements(this,elementtype,nb_elements)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_HybridElements), intent(inout) :: this
  type(atlas_mesh_ElementType) :: elementtype
  integer(c_size_t) :: nb_elements
  call atlas__mesh__HybridElements__add_elements(this%c_ptr(),elementtype%c_ptr(),nb_elements)
end subroutine

subroutine add_elements_with_nodes(this,elementtype,nb_elements,node_connectivity)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_HybridElements), intent(inout) :: this
  type(atlas_mesh_ElementType), intent(in) :: elementtype
  integer(c_size_t), intent(in) :: nb_elements
  integer(c_int), intent(in) :: node_connectivity(:)
  call atlas__mesh__HybridElements__add_elements_with_nodes(this%c_ptr(),&
    & elementtype%c_ptr(),nb_elements,node_connectivity,1)
end subroutine

! ===============================================================

subroutine atlas_mesh_ElementType__delete(this)
  use atlas_hybridelements_c_binding
  class(atlas_mesh_ElementType), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__mesh__ElementType__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_mesh_ElementType__copy(this,obj_in)
  class(atlas_mesh_ElementType), intent(inout) :: this
  class(atlas_RefCounted),   target, intent(in) :: obj_in
end subroutine


function atlas_mesh_Quadrilateral__constructor() result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_ElementType) :: this
  call this%reset_c_ptr( atlas__mesh__Quadrilateral__create() )
end function

function atlas_mesh_Triangle__constructor() result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_ElementType) :: this
  call this%reset_c_ptr( atlas__mesh__Triangle__create() )
end function

function atlas_mesh_Line__constructor() result(this)
  use atlas_hybridelements_c_binding
  type(atlas_mesh_ElementType) :: this
  call this%reset_c_ptr( atlas__mesh__Line__create() )
end function


end module atlas_mesh_HybridElements_module

