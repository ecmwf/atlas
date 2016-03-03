
module atlas_mesh_Nodes_module


use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr
use atlas_c_interop, only: c_str, c_to_f_string_cptr, atlas_free
use atlas_refcounted_module, only: atlas_refcounted
use atlas_Connectivity_module, only: atlas_Connectivity
use atlas_Field_module, only: atlas_Field
use atlas_Metadata_module, only: atlas_Metadata

implicit none

private :: c_size_t, c_int, c_ptr
private :: c_str, c_to_f_string_cptr, atlas_free
private :: atlas_refcounted
private :: atlas_Connectivity
private :: atlas_Field

public :: atlas_mesh_Nodes

private

!-----------------------------
! atlas_mesh_Nodes           !
!-----------------------------

TYPE, extends(atlas_refcounted) :: atlas_mesh_Nodes
contains
procedure, public :: size => atlas_mesh_Nodes__size
procedure, public :: resize
procedure, private :: add_field
procedure, private :: add_connectivity
generic, public :: add => &
    & add_field, &
    & add_connectivity
procedure, public :: remove_field
procedure, private :: field_by_idx_int
procedure, private :: field_by_idx_size_t
procedure, private :: field_by_name
generic, public :: field => &
    & field_by_idx_size_t, field_by_idx_int, &
    & field_by_name
procedure, public :: nb_fields
procedure, public :: has_field
procedure, public :: metadata
procedure, public :: str

procedure, public :: lonlat
procedure, public :: global_index
procedure, public :: remote_index
procedure, public :: partition
procedure, public :: ghost

procedure, public :: edge_connectivity
procedure, public :: cell_connectivity

procedure, public :: connectivity

procedure, public :: delete => atlas_mesh_Nodes__delete
procedure, public :: copy => atlas_mesh_Nodes__copy

end type

interface atlas_mesh_Nodes
  module procedure atlas_mesh_Nodes__cptr
  module procedure atlas_mesh_Nodes__constructor
end interface

!========================================================
contains
!========================================================

function atlas_mesh_Nodes__cptr(cptr) result(this)
  use atlas_nodes_c_binding
  type(atlas_mesh_Nodes) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
end function

function atlas_mesh_Nodes__constructor() result(this)
  use atlas_nodes_c_binding
  type(atlas_mesh_Nodes) :: this
  call this%reset_c_ptr( atlas__mesh__Nodes__create() )
end function

subroutine atlas_mesh_Nodes__delete(this)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__mesh__Nodes__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_mesh_Nodes__copy(this,obj_in)
  class(atlas_mesh_Nodes), intent(inout) :: this
  class(atlas_RefCounted),   target, intent(in) :: obj_in
end subroutine

function atlas_mesh_Nodes__size(this) result(val)
  use atlas_nodes_c_binding
  integer(c_size_t) :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  val = atlas__mesh__Nodes__size(this%c_ptr())
end function

function edge_connectivity(this) result(connectivity)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  type(atlas_Connectivity) :: connectivity
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__edge_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function cell_connectivity(this) result(connectivity)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  type(atlas_Connectivity) :: connectivity
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__cell_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function connectivity(this,name)
  use atlas_nodes_c_binding
  type(atlas_Connectivity) :: connectivity
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__connectivity(this%c_ptr(),c_str(name)) )
  call connectivity%return()
end function

subroutine add_connectivity(this,connectivity)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Connectivity), intent(in) :: connectivity
  call atlas__mesh__Nodes__add_connectivity(this%c_ptr(), connectivity%c_ptr())
end subroutine


subroutine add_field(this,field)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__mesh__Nodes__add_field(this%c_ptr(), field%c_ptr())
end subroutine

subroutine remove_field(this,name)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__mesh__Nodes__remove_field(this%c_ptr(),c_str(name))
end subroutine

function nb_fields(this) result(val)
  use atlas_nodes_c_binding
  integer(c_size_t) :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  val = atlas__mesh__Nodes__nb_fields(this%c_ptr())
end function

function has_field(this,name) result(val)
  use atlas_nodes_c_binding
  logical :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  if( atlas__mesh__Nodes__has_field(this%c_ptr(),c_str(name)) == 0 ) then
    val = .False.
  else
    val = .True.
  endif
end function

function field_by_name(this,name) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function field_by_idx_size_t(this,idx) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_size_t), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%c_ptr(),idx-1_c_size_t) )
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_int), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%c_ptr(),int(idx-1,c_size_t)) )
  call field%return()
end function

function lonlat(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__lonlat(this%c_ptr()) )
  call field%return()
end function

function global_index(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__global_index(this%c_ptr()) )
  call field%return()
end function

function remote_index(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__remote_index(this%c_ptr()) )
  call field%return()
end function

function partition(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__partition(this%c_ptr()) )
  call field%return()
end function

function ghost(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__ghost(this%c_ptr()) )
  call field%return()
end function

function metadata(this)
  use atlas_nodes_c_binding
  type(atlas_Metadata) :: metadata
  class(atlas_mesh_Nodes), intent(in) :: this
  call metadata%reset_c_ptr( atlas__mesh__Nodes__metadata(this%c_ptr()) )
end function

subroutine resize(this,size)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_size_t), intent(in) :: size
  call atlas__mesh__Nodes__resize(this%c_ptr(),size)
end subroutine

function str(this)
  use atlas_nodes_c_binding
  character(len=:), allocatable :: str
  class(atlas_mesh_Nodes), intent(in) :: this
  type(c_ptr) :: str_cptr
  integer(c_int) :: str_size
  call atlas__mesh__Nodes__str(this%c_ptr(),str_cptr,str_size)
  allocate(character(len=str_size) :: str )
  str = c_to_f_string_cptr(str_cptr)
  call atlas_free(str_cptr)
end function


! ----------------------------------------------------------------------------------------

end module atlas_mesh_Nodes_module

