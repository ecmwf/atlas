#include "atlas/atlas_f.h"

module atlas_mesh_Nodes_module

use fckit_owned_object_module, only: fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_mesh_Nodes

private

!-----------------------------
! atlas_mesh_Nodes           !
!-----------------------------

TYPE, extends(fckit_owned_object) :: atlas_mesh_Nodes
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

procedure, public :: xy
procedure, public :: lonlat
procedure, public :: global_index
procedure, public :: remote_index
procedure, public :: partition
procedure, public :: ghost

procedure, public :: edge_connectivity
procedure, public :: cell_connectivity

procedure, public :: connectivity

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_mesh_Nodes__final_auto
#endif
end type

interface atlas_mesh_Nodes
  module procedure atlas_mesh_Nodes__cptr
  module procedure atlas_mesh_Nodes__constructor
end interface

!========================================================
contains
!========================================================

function atlas_mesh_Nodes__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_nodes_c_binding
  type(atlas_mesh_Nodes) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_mesh_Nodes__constructor() result(this)
  use atlas_nodes_c_binding
  type(atlas_mesh_Nodes) :: this
  call this%reset_c_ptr( atlas__mesh__Nodes__create() )
  call this%return()
end function

function atlas_mesh_Nodes__size(this) result(val)
  use, intrinsic :: iso_c_binding, only: c_size_t
  use atlas_nodes_c_binding
  integer(c_size_t) :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  val = atlas__mesh__Nodes__size(this%c_ptr())
end function

function edge_connectivity(this) result(connectivity)
  use atlas_nodes_c_binding
  use atlas_Connectivity_module, only: atlas_Connectivity
  class(atlas_mesh_Nodes), intent(in) :: this
  type(atlas_Connectivity) :: connectivity
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__edge_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function cell_connectivity(this) result(connectivity)
  use atlas_nodes_c_binding
  use atlas_Connectivity_module, only: atlas_Connectivity
  class(atlas_mesh_Nodes), intent(in) :: this
  type(atlas_Connectivity) :: connectivity
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__cell_connectivity(this%c_ptr()) )
  call connectivity%return()
end function

function connectivity(this,name)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
  use atlas_Connectivity_module, only: atlas_Connectivity
  type(atlas_Connectivity) :: connectivity
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__connectivity(this%c_ptr(),c_str(name)) )
  call connectivity%return()
end function

subroutine add_connectivity(this,connectivity)
  use atlas_nodes_c_binding
  use atlas_Connectivity_module, only: atlas_Connectivity
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Connectivity), intent(in) :: connectivity
  call atlas__mesh__Nodes__add_connectivity(this%c_ptr(), connectivity%c_ptr())
end subroutine


subroutine add_field(this,field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__mesh__Nodes__add_field(this%c_ptr(), field%c_ptr())
end subroutine

subroutine remove_field(this,name)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__mesh__Nodes__remove_field(this%c_ptr(),c_str(name))
end subroutine

function nb_fields(this) result(val)
  use atlas_nodes_c_binding
  use, intrinsic :: iso_c_binding, only: c_size_t
  integer(c_size_t) :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  val = atlas__mesh__Nodes__nb_fields(this%c_ptr())
end function

function has_field(this,name) result(val)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
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
  use fckit_c_interop_module, only: c_str
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function field_by_idx_size_t(this,idx) result(field)
  use atlas_nodes_c_binding
  use, intrinsic :: iso_c_binding, only: c_size_t
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_size_t), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%c_ptr(),idx-1_c_size_t) )
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use atlas_nodes_c_binding
  use, intrinsic :: iso_c_binding, only: c_size_t, c_int
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_int), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%c_ptr(),int(idx-1,c_size_t)) )
  call field%return()
end function

function xy(this) result(field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__xy(this%c_ptr()) )
  call field%return()
end function

function lonlat(this) result(field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__lonlat(this%c_ptr()) )
  call field%return()
end function

function global_index(this) result(field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__global_index(this%c_ptr()) )
  call field%return()
end function

function remote_index(this) result(field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__remote_index(this%c_ptr()) )
  call field%return()
end function

function partition(this) result(field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__partition(this%c_ptr()) )
  call field%return()
end function

function ghost(this) result(field)
  use atlas_nodes_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__ghost(this%c_ptr()) )
  call field%return()
end function

function metadata(this)
  use atlas_nodes_c_binding
  use atlas_Metadata_module, only: atlas_Metadata
  type(atlas_Metadata) :: metadata
  class(atlas_mesh_Nodes), intent(in) :: this
  call metadata%reset_c_ptr( atlas__mesh__Nodes__metadata(this%c_ptr()) )
end function

subroutine resize(this,size)
  use atlas_nodes_c_binding
  use, intrinsic :: iso_c_binding, only: c_size_t
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_size_t), intent(in) :: size
  call atlas__mesh__Nodes__resize(this%c_ptr(),size)
end subroutine

function str(this)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_ptr_to_string, c_ptr_free

  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  character(len=:), allocatable :: str
  class(atlas_mesh_Nodes), intent(in) :: this
  type(c_ptr) :: str_cptr
  integer(c_int) :: str_size
  call atlas__mesh__Nodes__str(this%c_ptr(),str_cptr,str_size)
  allocate(character(len=str_size) :: str )
  str = c_ptr_to_string(str_cptr)
  call c_ptr_free(str_cptr)
end function

!-------------------------------------------------------------------------------

subroutine atlas_mesh_Nodes__final_auto(this)
  type(atlas_mesh_Nodes) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_mesh_Nodes__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_mesh_Nodes_module

