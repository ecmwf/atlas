
module atlas_ElementType_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, c_int
use fckit_c_interop_module, only : c_ptr_to_string
use fckit_refcounted_module, only : fckit_refcounted

implicit none

private :: c_ptr, c_size_t, c_int
private :: c_ptr_to_string
private :: fckit_refcounted

public :: atlas_ElementType
public :: atlas_Triangle
public :: atlas_Quadrilateral
public :: atlas_Line


private

!-----------------------------
! atlas_ElementType     !
!-----------------------------

type, extends(fckit_refcounted) :: atlas_ElementType
contains
! Public methods
  procedure, public :: copy     => atlas_ElementType__copy
  procedure, public :: delete   => atlas_ElementType__delete

  procedure, public :: nb_nodes
  procedure, public :: nb_edges
  procedure, public :: name
  procedure, public :: parametric

end type

interface atlas_ElementType
  module procedure atlas_ElementType__cptr
end interface

interface atlas_Triangle
  module procedure atlas_Triangle__constructor
end interface

interface atlas_Quadrilateral
  module procedure atlas_Quadrilateral__constructor
end interface

interface atlas_Line
  module procedure atlas_Line__constructor
end interface

!========================================================
contains
!========================================================

subroutine atlas_ElementType__delete(this)
  use atlas_elementtype_c_binding
  class(atlas_ElementType), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__mesh__ElementType__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_ElementType__copy(this,obj_in)
  class(atlas_ElementType), intent(inout) :: this
  class(fckit_refcounted),   target, intent(in) :: obj_in
end subroutine

function atlas_ElementType__cptr(cptr) result(this)
  use atlas_elementtype_c_binding
  type(atlas_ElementType) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
end function

function atlas_Quadrilateral__constructor() result(this)
  use atlas_elementtype_c_binding
  type(atlas_ElementType) :: this
  call this%reset_c_ptr( atlas__mesh__Quadrilateral__create() )
end function

function atlas_Triangle__constructor() result(this)
  use atlas_elementtype_c_binding
  type(atlas_ElementType) :: this
  call this%reset_c_ptr( atlas__mesh__Triangle__create() )
end function

function atlas_Line__constructor() result(this)
  use atlas_elementtype_c_binding
  type(atlas_ElementType) :: this
  call this%reset_c_ptr( atlas__mesh__Line__create() )
end function

function nb_nodes(this)
  use atlas_elementtype_c_binding
  integer(c_size_t) :: nb_nodes
  class(atlas_ElementType), intent(in) :: this
  nb_nodes = atlas__mesh__ElementType__nb_nodes(this%c_ptr())
end function

function nb_edges(this)
  use atlas_elementtype_c_binding
  integer(c_size_t) :: nb_edges
  class(atlas_ElementType), intent(in) :: this
  nb_edges = atlas__mesh__ElementType__nb_edges(this%c_ptr())
end function

function name(this)
  use atlas_elementtype_c_binding
  character(len=:), allocatable :: name
  class(atlas_ElementType) :: this
  type(c_ptr) :: name_c_str
  name_c_str = atlas__mesh__ElementType__name(this%c_ptr())
  name = c_ptr_to_string(name_c_str)
end function

function parametric(this)
  use atlas_elementtype_c_binding
  logical :: parametric
  class(atlas_ElementType), intent(in) :: this
  integer(c_int) :: parametric_int
  parametric_int = atlas__mesh__ElementType__parametric(this%c_ptr())
  if( parametric_int == 0 ) then
    parametric = .False.
  else
    parametric = .True.
  endif
end function


end module atlas_ElementType_module

