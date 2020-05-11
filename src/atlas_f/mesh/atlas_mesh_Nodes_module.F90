! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_mesh_Nodes_module

use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_long
use fckit_owned_object_module, only: fckit_owned_object
use atlas_Metadata_module, only: atlas_Metadata
use atlas_Field_module, only: atlas_Field
use atlas_Connectivity_module, only: atlas_Connectivity
use atlas_kinds_module, only : ATLAS_KIND_IDX

implicit none

private :: fckit_owned_object
private :: c_ptr, c_int
private :: atlas_Metadata
private :: atlas_Field
private :: atlas_Connectivity
public :: atlas_mesh_Nodes

private

!-----------------------------
! atlas_mesh_Nodes           !
!-----------------------------

TYPE, extends(fckit_owned_object) :: atlas_mesh_Nodes
contains
procedure, public :: size => atlas_mesh_Nodes__size
procedure, private :: resize_int
procedure, private :: resize_long
generic, public :: resize => resize_int, resize_long
procedure, private :: add_field
procedure, private :: add_connectivity
generic, public :: add => &
    & add_field, &
    & add_connectivity
procedure, public :: remove_field
procedure, private :: field_by_idx_int
procedure, private :: field_by_idx_long
procedure, private :: field_by_name
generic, public :: field => &
    & field_by_idx_long, field_by_idx_int, &
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
  use atlas_nodes_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  val = atlas__mesh__Nodes__size(this%CPTR_PGIBUG_A)
end function

function edge_connectivity(this) result(connectivity)
  use atlas_nodes_c_binding
  use atlas_Connectivity_module, only: atlas_Connectivity
  class(atlas_mesh_Nodes), intent(in) :: this
  type(atlas_Connectivity) :: connectivity
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__edge_connectivity(this%CPTR_PGIBUG_A) )
  call connectivity%return()
end function

function cell_connectivity(this) result(connectivity)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  type(atlas_Connectivity) :: connectivity
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__cell_connectivity(this%CPTR_PGIBUG_A) )
  call connectivity%return()
end function

function connectivity(this,name)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
  type(atlas_Connectivity) :: connectivity
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  connectivity = atlas_Connectivity( &
      atlas__mesh__Nodes__connectivity(this%CPTR_PGIBUG_A,c_str(name)) )
  call connectivity%return()
end function

subroutine add_connectivity(this,connectivity)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Connectivity), intent(in) :: connectivity
  call atlas__mesh__Nodes__add_connectivity(this%CPTR_PGIBUG_A, connectivity%CPTR_PGIBUG_A)
end subroutine


subroutine add_field(this,field)
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(inout) :: this
  type(atlas_Field), intent(in) :: field
  call atlas__mesh__Nodes__add_field(this%CPTR_PGIBUG_A, field%CPTR_PGIBUG_A)
end subroutine

subroutine remove_field(this,name)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__mesh__Nodes__remove_field(this%CPTR_PGIBUG_A,c_str(name))
end subroutine

function nb_fields(this) result(val)
  use atlas_nodes_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  val = atlas__mesh__Nodes__nb_fields(this%CPTR_PGIBUG_A)
end function

function has_field(this,name) result(val)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
  logical :: val
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  if( atlas__mesh__Nodes__has_field(this%CPTR_PGIBUG_A,c_str(name)) == 0 ) then
    val = .False.
  else
    val = .True.
  endif
end function

function field_by_name(this,name) result(field)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_str
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__mesh__Nodes__field_by_name(this%CPTR_PGIBUG_A,c_str(name)) )
  call field%return()
end function

function field_by_idx_long(this,idx) result(field)
  use atlas_nodes_c_binding
  use, intrinsic :: iso_c_binding, only: c_long
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_long), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%CPTR_PGIBUG_A,int(idx-1_c_long,ATLAS_KIND_IDX) ) )
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use atlas_nodes_c_binding
  use, intrinsic :: iso_c_binding, only: c_int
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_int), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Nodes__field_by_idx(this%CPTR_PGIBUG_A,int(idx-1_c_long,ATLAS_KIND_IDX) ) )
  call field%return()
end function

function xy(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__xy(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function lonlat(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__lonlat(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function global_index(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__global_index(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function remote_index(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__remote_index(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function partition(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__partition(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function ghost(this) result(field)
  use atlas_nodes_c_binding
  type(atlas_Field) :: field
  class(atlas_mesh_Nodes), intent(in) :: this
  field = atlas_Field( atlas__mesh__Nodes__ghost(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function metadata(this)
  use atlas_nodes_c_binding
  type(atlas_Metadata) :: metadata
  class(atlas_mesh_Nodes), intent(in) :: this
  call metadata%reset_c_ptr( atlas__mesh__Nodes__metadata(this%CPTR_PGIBUG_A) )
end function

subroutine resize_int(this,size)
  use, intrinsic :: iso_c_binding
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_int), intent(in) :: size
  call atlas__mesh__Nodes__resize(this%CPTR_PGIBUG_A,int(size,ATLAS_KIND_IDX))
end subroutine

subroutine resize_long(this,size)
  use, intrinsic :: iso_c_binding
  use atlas_nodes_c_binding
  class(atlas_mesh_Nodes), intent(in) :: this
  integer(c_long), intent(in) :: size
  call atlas__mesh__Nodes__resize(this%CPTR_PGIBUG_A,int(size,ATLAS_KIND_IDX))
end subroutine

function str(this)
  use atlas_nodes_c_binding
  use fckit_c_interop_module, only: c_ptr_to_string, c_ptr_free
  character(len=:), allocatable :: str
  class(atlas_mesh_Nodes), intent(in) :: this
  type(c_ptr) :: str_cptr
  integer(c_int) :: str_size
  call atlas__mesh__Nodes__str(this%CPTR_PGIBUG_A,str_cptr,str_size)
  allocate(character(len=str_size) :: str )
  str = c_ptr_to_string(str_cptr)
  call c_ptr_free(str_cptr)
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_mesh_Nodes__final_auto(this)
  type(atlas_mesh_Nodes), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_mesh_Nodes__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! ----------------------------------------------------------------------------------------

end module atlas_mesh_Nodes_module

