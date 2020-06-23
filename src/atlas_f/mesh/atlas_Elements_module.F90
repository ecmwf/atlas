! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Elements_module

use fckit_owned_object_module, only: fckit_owned_object
use atlas_kinds_module, only : ATLAS_KIND_IDX

implicit none

private :: fckit_owned_object

public :: atlas_Elements

private

!-----------------------------
! atlas_Elements        !
!-----------------------------

type, extends(fckit_owned_object) :: atlas_Elements
contains
! Public methods
  procedure, public :: size     => atlas_Elements__size
  procedure, public :: begin => atlas_Elements__begin
  procedure, public :: end   => atlas_Elements__end

  procedure, public ::  node_connectivity
  procedure, public ::  edge_connectivity
  procedure, public ::  cell_connectivity

  generic, public :: add   => add_elements_long, add_elements_int
  generic, public :: field => field_by_idx_long, field_by_idx_int, field_by_name

  procedure, public :: element_type

  procedure, public :: nb_fields
  procedure, public :: has_field
  procedure, public :: global_index
  procedure, public :: remote_index
  procedure, public :: partition
  procedure, public :: halo

! Private methods
  procedure, private :: field_by_idx_int
  procedure, private :: field_by_idx_long
  procedure, private :: field_by_name
  procedure, private :: add_elements_long
  procedure, private :: add_elements_int

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Elements__final_auto
#endif
end type

interface atlas_Elements
  module procedure atlas_Elements__cptr
end interface


!========================================================
contains
!========================================================

function atlas_Elements__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_Elements) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_Elements__size(this) result(val)
  use atlas_elements_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__size(this%CPTR_PGIBUG_A)
end function

function node_connectivity(this) result(connectivity)
  use atlas_elements_c_binding
  use atlas_Connectivity_module, only: atlas_BlockConnectivity
  class(atlas_Elements), intent(in) :: this
  type(atlas_BlockConnectivity) :: connectivity
  connectivity = atlas_BlockConnectivity( &
      atlas__mesh__Elements__node_connectivity(this%CPTR_PGIBUG_A) )
  !call connectivity%return()
end function

function edge_connectivity(this) result(connectivity)
  use atlas_elements_c_binding
  use atlas_Connectivity_module, only: atlas_BlockConnectivity
  class(atlas_Elements), intent(in) :: this
  type(atlas_BlockConnectivity) :: connectivity
  connectivity = atlas_BlockConnectivity( &
      atlas__mesh__Elements__edge_connectivity(this%CPTR_PGIBUG_A) )
  !call connectivity%return()
end function

function cell_connectivity(this) result(connectivity)
  use atlas_elements_c_binding
  use atlas_Connectivity_module, only: atlas_BlockConnectivity
  class(atlas_Elements), intent(in) :: this
  type(atlas_BlockConnectivity) :: connectivity
  connectivity = atlas_BlockConnectivity( &
      atlas__mesh__Elements__cell_connectivity(this%CPTR_PGIBUG_A) )
  !call connectivity%return()
end function

function element_type(this) result(etype)
  use atlas_elements_c_binding
  use atlas_ElementType_module, only: atlas_ElementType
  class(atlas_Elements), intent(in) :: this
  type(atlas_ElementType) :: etype
  etype = atlas_ElementType( &
      atlas__mesh__Elements__element_type(this%CPTR_PGIBUG_A) )
  call etype%return()
end function

subroutine add_elements_long(this,nb_elements)
  use, intrinsic :: iso_c_binding, only: c_long, c_int
  use atlas_elements_c_binding
  class(atlas_Elements), intent(inout) :: this
  integer(c_long) :: nb_elements
  call atlas__mesh__Elements__add(this%CPTR_PGIBUG_A,int(nb_elements,ATLAS_KIND_IDX))
end subroutine

subroutine add_elements_int(this,nb_elements)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_elements_c_binding
  class(atlas_Elements), intent(inout) :: this
  integer(c_int) :: nb_elements
  call atlas__mesh__Elements__add(this%CPTR_PGIBUG_A,int(nb_elements,ATLAS_KIND_IDX))
end subroutine

function nb_fields(this) result(val)
  use atlas_elements_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__nb_fields(this%CPTR_PGIBUG_A)
end function

function has_field(this,name) result(val)
  use fckit_c_interop_module, only: c_str
  use atlas_elements_c_binding
  logical :: val
  class(atlas_Elements), intent(in) :: this
  character(len=*), intent(in) :: name
  if( atlas__mesh__Elements__has_field(this%CPTR_PGIBUG_A,c_str(name)) == 0 ) then
    val = .False.
  else
    val = .True.
  endif
end function

function field_by_name(this,name) result(field)
  use fckit_c_interop_module, only: c_str
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__mesh__Elements__field_by_name(this%CPTR_PGIBUG_A,c_str(name)) )
  call field%return()
end function

function field_by_idx_long(this,idx) result(field)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  integer(c_long), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Elements__field_by_idx(this%CPTR_PGIBUG_A,int(idx-1_c_long,ATLAS_KIND_IDX) ) )
  call field%return()
end function

function field_by_idx_int(this,idx) result(field)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  integer(c_int), intent(in) :: idx
  field = atlas_Field( atlas__mesh__Elements__field_by_idx(this%CPTR_PGIBUG_A,int(idx-1_c_int,ATLAS_KIND_IDX) ) )
  call field%return()
end function

function global_index(this) result(field)
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__global_index(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function remote_index(this) result(field)
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__remote_index(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function partition(this) result(field)
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__partition(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function halo(this) result(field)
  use atlas_elements_c_binding
  use atlas_Field_module, only: atlas_Field
  type(atlas_Field) :: field
  class(atlas_Elements), intent(in) :: this
  field = atlas_Field( atlas__mesh__Elements__halo(this%CPTR_PGIBUG_A) )
  call field%return()
end function

function atlas_Elements__begin(this) result(val)
  use atlas_elements_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__begin(this%CPTR_PGIBUG_A) + 1
end function

function atlas_Elements__end(this) result(val)
  use atlas_elements_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Elements), intent(in) :: this
  val = atlas__mesh__Elements__end(this%CPTR_PGIBUG_A)
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Elements__final_auto(this)
  type(atlas_Elements), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Elements__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_Elements_module
