
module atlas_Elements_module

!use, intrinsic :: iso_c_binding, only : c_funptr, c_ptr, c_loc, c_f_pointer, c_f_procpointer, c_funloc, c_int, c_size_t
use fckit_refcounted_module, only: fckit_refcounted
use atlas_Connectivity_module, only: atlas_BlockConnectivity
use atlas_ElementType_module, only: atlas_ElementType
use atlas_Field_module, only: atlas_Field
use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr
use fckit_c_interop_module, only: c_str

implicit none

private :: c_size_t, c_int, c_ptr
private :: c_str
private :: fckit_refcounted
private :: atlas_BlockConnectivity
private :: atlas_Field
private :: atlas_ElementType

public :: atlas_Elements

private

!-----------------------------
! atlas_Elements        !
!-----------------------------

type, extends(fckit_refcounted) :: atlas_Elements
contains
! Public methods
  procedure, public :: delete   => atlas_Elements__delete

  procedure, public :: size     => atlas_Elements__size
  procedure, public :: begin => atlas_Elements__begin
  procedure, public :: end   => atlas_Elements__end

  procedure, public ::  node_connectivity
  procedure, public ::  edge_connectivity
  procedure, public ::  cell_connectivity

  generic, public :: add   => add_elements_size_t, add_elements_int
  generic, public :: field => field_by_idx_size_t, field_by_idx_int, field_by_name

  procedure, public :: element_type

  procedure, public :: nb_fields
  procedure, public :: has_field
  procedure, public :: global_index
  procedure, public :: remote_index
  procedure, public :: partition
  procedure, public :: halo

! Private methods
  procedure, private :: field_by_idx_int
  procedure, private :: field_by_idx_size_t
  procedure, private :: field_by_name
  procedure, private :: add_elements_size_t
  procedure, private :: add_elements_int

end type

interface atlas_Elements
  module procedure atlas_Elements__cptr
end interface


!========================================================
contains
!========================================================

function atlas_Elements__cptr(cptr) result(elements)
  type(atlas_Elements) :: elements
  type(c_ptr), intent(in) :: cptr
  call elements%reset_c_ptr( cptr )
end function

subroutine atlas_Elements__delete(this)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__mesh__Elements__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

function atlas_Elements__size(this) result(val)
  use atlas_elements_c_binding
  integer(c_size_t) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__size(this%c_ptr())
end function

function node_connectivity(this) result(connectivity)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(in) :: this
  type(atlas_BlockConnectivity) :: connectivity
  connectivity = atlas_BlockConnectivity( &
      atlas__mesh__Elements__node_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function edge_connectivity(this) result(connectivity)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(in) :: this
  type(atlas_BlockConnectivity) :: connectivity
  connectivity = atlas_BlockConnectivity( &
      atlas__mesh__Elements__edge_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function cell_connectivity(this) result(connectivity)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(in) :: this
  type(atlas_BlockConnectivity) :: connectivity
  connectivity = atlas_BlockConnectivity( &
      atlas__mesh__Elements__cell_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function element_type(this) result(etype)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(in) :: this
  type(atlas_ElementType) :: etype
  etype = atlas_ElementType( &
      atlas__mesh__Elements__element_type(this%c_ptr()) )
  call etype%return()
end function

subroutine add_elements_size_t(this,nb_elements)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(inout) :: this
  integer(c_size_t) :: nb_elements
  call atlas__mesh__Elements__add(this%c_ptr(),nb_elements)
end subroutine

subroutine add_elements_int(this,nb_elements)
  use atlas_elements_c_binding
  class(atlas_Elements), intent(inout) :: this
  integer(c_int) :: nb_elements
  call atlas__mesh__Elements__add(this%c_ptr(),int(nb_elements,c_size_t))
end subroutine

function nb_fields(this) result(val)
  use atlas_elements_c_binding
  integer(c_size_t) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__nb_fields(this%c_ptr())
end function

function has_field(this,name) result(val)
  use atlas_elements_c_binding
  logical :: val
  class(atlas_Elements), intent(in) :: this
  character(len=*), intent(in) :: name
  if( atlas__mesh__Elements__has_field(this%c_ptr(),c_str(name)) == 0 ) then
    val = .False.
  else
    val = .True.
  endif
end function

function field_by_name(this,name) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__mesh__Elements__field_by_name(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function field_by_idx_size_t(this,idx) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  integer(c_size_t), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Elements__field_by_idx(this%c_ptr(),idx-1_c_size_t) )
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  integer(c_int), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Elements__field_by_idx(this%c_ptr(),int(idx-1,c_size_t)) )
  call field%return()
end function

function global_index(this) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__global_index(this%c_ptr()) )
  call field%return()
end function

function remote_index(this) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__remote_index(this%c_ptr()) )
  call field%return()
end function

function partition(this) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__partition(this%c_ptr()) )
  call field%return()
end function

function halo(this) result(field)
  use atlas_elements_c_binding
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__halo(this%c_ptr()) )
  call field%return()
end function

function atlas_Elements__begin(this) result(val)
  use atlas_elements_c_binding
  integer(c_size_t) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__begin(this%c_ptr()) + 1
end function

function atlas_Elements__end(this) result(val)
  use atlas_elements_c_binding
  integer(c_size_t) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__end(this%c_ptr())
end function

end module atlas_Elements_module
