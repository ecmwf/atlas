! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_connectivity_module

use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr
use fckit_owned_object_module, only : fckit_owned_object
use fckit_object_module, only : fckit_object
use atlas_kinds_module, only : ATLAS_KIND_IDX
use atlas_allocate_module, only : atlas_allocate_managedmem, atlas_deallocate_managedmem
use fckit_exception_module, only : fckit_exception
implicit none

private :: c_ptr
private :: c_int
private :: c_null_ptr
private :: fckit_owned_object
private :: fckit_object
private :: atlas_allocate_managedmem
private :: atlas_deallocate_managedmem

public :: atlas_Connectivity
public :: atlas_MultiBlockConnectivity
public :: atlas_BlockConnectivity


private

!-----------------------------
! atlas_Connectivity         !
!-----------------------------

type, extends(fckit_owned_object) :: atlas_Connectivity

! Public members
  type( atlas_ConnectivityAccess ), public, pointer :: access => null()

contains

  procedure, public :: assignment_operator_hook
  procedure, public :: set_access

! Public methods
  procedure, public  :: name => atlas_Connectivity__name
  procedure, private :: value_args_int     => atlas_Connectivity__value_args_int
  procedure, private :: value_args_long    => atlas_Connectivity__value_args_long
  generic, public :: value => value_args_int, value_args_long
  procedure, public :: rows     => atlas_Connectivity__rows
  procedure, public :: cols     => atlas_Connectivity__cols
  procedure, public :: maxcols  => atlas_Connectivity__maxcols
  procedure, public :: mincols  => atlas_Connectivity__mincols
  procedure, public :: padded_data     => atlas_Connectivity__padded_data
  procedure, private :: atlas_Connectivity__data
  procedure, private :: atlas_Connectivity__data_int
  procedure, private :: atlas_Connectivity__data_long
  generic, public :: data     => atlas_Connectivity__data, atlas_Connectivity__data_int, atlas_Connectivity__data_long
  procedure, private :: atlas_Connectivity__row_int
  procedure, private :: atlas_Connectivity__row_long
  generic, public :: row      => atlas_Connectivity__row_int, atlas_Connectivity__row_long
  procedure, public :: missing_value  => atlas_Connectivity__missing_value
  procedure, private :: add_values_args_idx => atlas_Connectivity__add_values_args_idx
#if ATLAS_BITS_LOCAL != 32
  procedure, private :: add_values_args_int32 => atlas_Connectivity__add_values_args_int32
#define _add_values_args_int32  , add_values_args_int32
#else
#define _add_values_args_int32
#endif
  procedure, private :: add_values_args_long => atlas_Connectivity__add_values_args_long
  procedure, private :: add_missing_args_int => atlas_Connectivity__add_missing_args_int
  procedure, private :: add_missing_args_long => atlas_Connectivity__add_missing_args_long
  generic, public :: add => add_values_args_idx, add_values_args_long, add_missing_args_int, &
    & add_missing_args_long _add_values_args_int32
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Connectivity__final_auto
#endif

end type

!---------------------------------------
! atlas_MultiBlockConnectivity         !
!---------------------------------------

type, extends(atlas_Connectivity) :: atlas_MultiBlockConnectivity
contains
! Public methods
  procedure, public :: blocks   => atlas_MultiBlockConnectivity__blocks
  procedure, public :: block    => atlas_MultiBlockConnectivity__block

! PGI compiler bug won't accept "assignment_operator_hook" from atlas_Connectivity parent class... grrr
  procedure, public :: assignment_operator_hook => atlas_MBC__assignment_operator_hook

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_MultiBlockConnectivity__final_auto
#endif

end type

!----------------------------!
! atlas_BlockConnectivity    !
!----------------------------!

type, extends(fckit_object) :: atlas_BlockConnectivity
contains
  procedure, public :: rows     => atlas_BlockConnectivity__rows
  procedure, public :: cols     => atlas_BlockConnectivity__cols
  procedure, public :: data     => atlas_BlockConnectivity__data
  procedure, public :: missing_value  => atlas_BlockConnectivity__missing_value
end type

interface atlas_Connectivity
  module procedure Connectivity_cptr
  module procedure Connectivity_constructor
end interface

interface atlas_MultiBlockConnectivity
  module procedure MultiBlockConnectivity_cptr
  module procedure MultiBlockConnectivity_constructor
end interface

interface atlas_BlockConnectivity
  module procedure BlockConnectivity_cptr
end interface


!-------------------------------
! Helper types                 !
!-------------------------------

type :: atlas_ConnectivityAccessRow
  integer(ATLAS_KIND_IDX), pointer :: col(:)
contains
end type

type :: atlas_ConnectivityAccess
  integer(ATLAS_KIND_IDX), private, pointer :: values_(:) => null()
  integer(ATLAS_KIND_IDX), private, pointer :: displs_(:) => null()
  integer(ATLAS_KIND_IDX), public,  pointer :: cols(:)    => null()
  type(atlas_ConnectivityAccessRow), public, pointer :: row(:) => null()
  integer(ATLAS_KIND_IDX), private :: rows_
  integer(ATLAS_KIND_IDX), private, pointer :: padded_(:,:) => null()
  integer(ATLAS_KIND_IDX), private :: maxcols_, mincols_
  integer(ATLAS_KIND_IDX), private :: missing_value_
  type(c_ptr), private :: connectivity_ptr = c_null_ptr
contains
  procedure, public :: rows  => access_rows
  procedure, public :: value => access_value
end type

! ----------------------------------------------------------------------------
interface
    subroutine atlas__connectivity__register_ctxt(This,ptr) bind(c,name="atlas__connectivity__register_ctxt")
        use, intrinsic :: iso_c_binding, only: c_ptr
        type(c_ptr),    value :: This
        type(c_ptr),    value :: ptr
    end subroutine

    subroutine atlas__connectivity__register_update(This,funptr) bind(c,name="atlas__connectivity__register_update")
        use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr
        type(c_ptr),    value :: This
        type(c_funptr), value :: funptr
    end subroutine

    subroutine atlas__connectivity__register_delete(This,funptr) bind(c,name="atlas__connectivity__register_delete")
        use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr
        type(c_ptr),    value :: This
        type(c_funptr), value :: funptr
    end subroutine

    function atlas__connectivity__ctxt(This,ptr) bind(c,name="atlas__connectivity__ctxt")
        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        integer(c_int) :: atlas__connectivity__ctxt
        type(c_ptr),    value :: This
        type(c_ptr) :: ptr
    end function

end interface

contains

subroutine assignment_operator_hook(this,other)
  class(atlas_Connectivity) :: this
  class(fckit_owned_object) :: other
  call this%set_access()
  FCKIT_SUPPRESS_UNUSED(other)
end subroutine

! Following routine is exact copy of "assignment_operator_hook" above, because of bug in PGI compiler (17.7)
! Without it, wrongly the "fckit_owned_object::assignment_operator_hook" is used instead of
! "atlas_Connectivity::assignment_operator_hook".
subroutine atlas_MBC__assignment_operator_hook(this,other)
  class(atlas_MultiBlockConnectivity) :: this
  class(fckit_owned_object) :: other
  call this%set_access()
  FCKIT_SUPPRESS_UNUSED(other)
end subroutine

function Connectivity_cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_connectivity_c_binding
  type(atlas_Connectivity) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
  call this%set_access()
  call this%return()
end function

function Connectivity_constructor(name) result(this)
  use atlas_connectivity_c_binding
  use fckit_c_interop_module, only : c_str
  type(atlas_Connectivity) :: this
  character(len=*), intent(in), optional :: name
  call this%reset_c_ptr( atlas__Connectivity__create() )
  call this%set_access()
  if( present(name) ) then
    call atlas__Connectivity__rename(this%CPTR_PGIBUG_A,c_str(name))
  endif
  call this%return()
end function


function atlas_Connectivity__name(this) result(name)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_connectivity_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string
  class(atlas_Connectivity), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__Connectivity__name(this%CPTR_PGIBUG_A)
  name = c_ptr_to_string(name_c_str)
end function

pure function access_value(this,c,r) result(val)
  use, intrinsic :: iso_c_binding, only : c_int
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_ConnectivityAccess), intent(in) :: this
  integer(ATLAS_KIND_IDX), intent(in) :: r,c
  val = this%values_(c+this%displs_(r))
end function

pure function access_rows(this) result(val)
  use atlas_kinds_module
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_ConnectivityAccess), intent(in) :: this
  val = this%rows_
end function

pure function atlas_Connectivity__value_args_long(this,c,r) result(val)
  use, intrinsic :: iso_c_binding, only : c_long
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_long), intent(in) :: r,c
  val = this%access%values_(c+this%access%displs_(r))
end function

pure function atlas_Connectivity__value_args_int(this,c,r) result(val)
  use, intrinsic :: iso_c_binding, only : c_int
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: r,c
  val = this%access%values_(c+this%access%displs_(r))
end function

pure function atlas_Connectivity__rows(this) result(val)
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%rows_
end function

function atlas_Connectivity__missing_value(this) result(val)
  use atlas_connectivity_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = atlas__Connectivity__missing_value(this%CPTR_PGIBUG_A)
end function

pure function atlas_Connectivity__cols(this,r) result(val)
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: r
  val = this%access%cols(r)
end function

pure function atlas_Connectivity__mincols(this) result(val)
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%mincols_
end function

pure function atlas_Connectivity__maxcols(this) result(val)
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%maxcols_
end function

subroutine atlas_Connectivity__padded_data(this, padded, cols)
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_Connectivity), intent(inout) :: this
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: padded(:,:)
  integer(ATLAS_KIND_IDX), pointer, intent(inout), optional :: cols(:)
  if( .not. associated(this%access%padded_) ) call update_padded(this%access)
  padded => this%access%padded_
  if( present(cols) ) cols => this%access%cols
end subroutine

function c_loc_idx(x)
  use, intrinsic :: iso_c_binding
  integer(ATLAS_KIND_IDX), target :: x
  type(c_ptr) :: c_loc_idx
  c_loc_idx = c_loc(x)
end function

subroutine atlas_Connectivity__data(this, data)
  use, intrinsic :: iso_c_binding, only : c_f_pointer, c_loc, c_int
  class(atlas_Connectivity), intent(in) :: this
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: data(:,:)
  integer(ATLAS_KIND_IDX) :: maxcols

  maxcols = this%maxcols()
  if( maxcols == this%mincols() ) then
    if( size(this%access%values_) > 0 ) then
      call c_f_pointer (c_loc_idx(this%access%values_(1)), data, &
          & (/maxcols,this%access%rows_/))
    endif
  else
    call fckit_exception%abort("ERROR: Cannot point connectivity pointer data(:,:) to irregular table",&
      & "atlas_Connectivity_module.F90", __LINE__ )
  endif
end subroutine

subroutine atlas_Connectivity__data_int(this, data, ncols)
  use, intrinsic :: iso_c_binding, only : c_f_pointer, c_loc, c_int
  class(atlas_Connectivity), intent(in) :: this
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: data(:,:)
  integer(c_int), intent(out) :: ncols
  integer(ATLAS_KIND_IDX) :: maxcols

  maxcols = this%maxcols()
  if( maxcols == this%mincols() ) then
    if( size(this%access%values_) > 0 ) then
      call c_f_pointer (c_loc_idx(this%access%values_(1)), data, &
          & (/maxcols,this%access%rows_/))
      ncols = maxcols
    endif
  else
    call fckit_exception%abort("ERROR: Cannot point connectivity pointer data(:,:) to irregular table",&
      & "atlas_Connectivity_module.F90", __LINE__ )
  endif
end subroutine

subroutine atlas_Connectivity__data_long(this, data, ncols)
  use, intrinsic :: iso_c_binding, only : c_f_pointer, c_loc, c_long
  class(atlas_Connectivity), intent(in) :: this
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: data(:,:)
  integer(c_long), intent(out) :: ncols
  integer(ATLAS_KIND_IDX) :: maxcols

  maxcols = this%maxcols()
  if( maxcols == this%mincols() ) then
    if( size(this%access%values_) > 0 ) then
      call c_f_pointer (c_loc_idx(this%access%values_(1)), data, &
          & (/maxcols,this%access%rows_/))
      ncols = maxcols
    endif
  else
    call fckit_exception%abort("ERROR: Cannot point connectivity pointer data(:,:) to irregular table",&
      & "atlas_Connectivity_module.F90", __LINE__ )
  endif
end subroutine

subroutine atlas_Connectivity__row_int(this, row_idx, row, cols)
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: row_idx
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: row(:)
  integer(c_int), intent(out) ::  cols
  row  => this%access%values_(this%access%displs_(row_idx)+1:this%access%displs_(row_idx+1)+1)
  cols =  this%access%cols(row_idx)
end subroutine

subroutine atlas_Connectivity__row_long(this, row_idx, row, cols)
  use, intrinsic :: iso_c_binding, only : c_long
  class(atlas_Connectivity), intent(in) :: this
  integer(c_long), intent(in) :: row_idx
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: row(:)
  integer(c_long), intent(out) ::  cols
  row  => this%access%values_(this%access%displs_(row_idx)+1:this%access%displs_(row_idx+1)+1)
  cols =  this%access%cols(row_idx)
end subroutine

subroutine atlas_Connectivity__add_values_args_long(this,rows,cols,values)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_long
  class(atlas_Connectivity), intent(in) :: this
  integer(c_long), intent(in) :: rows
  integer(c_long), intent(in) :: cols
  integer(ATLAS_KIND_IDX) :: values(:)
  call atlas__connectivity__add_values(this%CPTR_PGIBUG_A,int(rows,ATLAS_KIND_IDX),int(cols,ATLAS_KIND_IDX),values)
end subroutine

subroutine atlas_Connectivity__add_values_args_idx(this,rows,cols,values)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: rows
  integer(c_int), intent(in) :: cols
  integer(ATLAS_KIND_IDX), intent(in) :: values(:)
  call atlas__connectivity__add_values(this%CPTR_PGIBUG_A,int(rows,ATLAS_KIND_IDX),int(cols,ATLAS_KIND_IDX),values)
end subroutine

#if ATLAS_BITS_LOCAL != 32
subroutine atlas_Connectivity__add_values_args_int32(this,rows,cols,values)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: rows
  integer(c_int), intent(in) :: cols
  integer(c_int), intent(in) :: values(:)
  integer(ATLAS_KIND_IDX) :: idx_values(rows*cols)
  idx_values(:) = values(:)
  call atlas__connectivity__add_values(this%CPTR_PGIBUG_A,int(rows,ATLAS_KIND_IDX),int(cols,ATLAS_KIND_IDX),idx_values)
end subroutine
#endif


subroutine atlas_Connectivity__add_missing_args_long(this,rows,cols)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_long
  class(atlas_Connectivity), intent(in) :: this
  integer(c_long) :: rows
  integer(c_long) :: cols
  call atlas__connectivity__add_missing(this%CPTR_PGIBUG_A,int(rows,ATLAS_KIND_IDX),int(cols,ATLAS_KIND_IDX))
end subroutine

subroutine atlas_Connectivity__add_missing_args_int(this,rows,cols)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_int
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int) :: rows
  integer(c_int) :: cols
  call atlas__connectivity__add_missing(this%CPTR_PGIBUG_A,int(rows,ATLAS_KIND_IDX),int(cols,ATLAS_KIND_IDX))
end subroutine

!========================================================

function MultiBlockConnectivity_cptr(cptr) result(this)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_MultiBlockConnectivity) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
  call this%set_access()
  call this%return()
end function

function MultiBlockConnectivity_constructor(name) result(this)
  use atlas_connectivity_c_binding
  use fckit_c_interop_module
  type(atlas_MultiBlockConnectivity) :: this
  character(len=*), intent(in), optional :: name
  call this%reset_c_ptr( atlas__MultiBlockConnectivity__create() )
  call this%set_access()
  if( present(name) ) then
    call atlas__Connectivity__rename(this%CPTR_PGIBUG_A,c_str(name))
  endif
  call this%return()
end function

function atlas_MultiBlockConnectivity__blocks(this) result(val)
  use atlas_connectivity_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_MultiBlockConnectivity), intent(in) :: this
  val = atlas__MultiBlockConnectivity__blocks(this%CPTR_PGIBUG_A)
end function

function atlas_MultiBlockConnectivity__block(this,block_idx) result(block)
  use atlas_connectivity_c_binding
  type(atlas_BlockConnectivity) :: block
  class(atlas_MultiBlockConnectivity), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: block_idx
  call block%reset_c_ptr( atlas__MultiBlockConnectivity__block( &
    this%CPTR_PGIBUG_A,int(block_idx-1_ATLAS_KIND_IDX,ATLAS_KIND_IDX) ) )
end function

!========================================================

function BlockConnectivity_cptr(cptr) result(this)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_BlockConnectivity) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
end function

subroutine atlas_BlockConnectivity__data(this,data)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
  class(atlas_BlockConnectivity), intent(in) :: this
  integer(ATLAS_KIND_IDX), pointer, intent(inout) :: data(:,:)
  type(c_ptr) :: data_cptr
  integer(ATLAS_KIND_IDX) :: rows
  integer(ATLAS_KIND_IDX) :: cols
  call atlas__BlockConnectivity__data(this%CPTR_PGIBUG_A,data_cptr,rows,cols)
  call c_f_pointer (data_cptr, data, [cols,rows])
end subroutine

function atlas_BlockConnectivity__rows(this) result(val)
  use atlas_connectivity_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_BlockConnectivity), intent(in) :: this
  val = atlas__BlockConnectivity__rows(this%CPTR_PGIBUG_A)
end function

function atlas_BlockConnectivity__cols(this) result(val)
  use atlas_connectivity_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_BlockConnectivity), intent(in) :: this
  val = atlas__BlockConnectivity__cols(this%CPTR_PGIBUG_A)
end function

function atlas_BlockConnectivity__missing_value(this) result(val)
  use atlas_connectivity_c_binding
  integer(ATLAS_KIND_IDX) :: val
  class(atlas_BlockConnectivity), intent(in) :: this
  val = atlas__BlockConnectivity__missing_value(this%CPTR_PGIBUG_A) + 1
end function

!========================================================

subroutine set_access(this)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_f_pointer, c_funloc, c_loc
  class(atlas_Connectivity), intent(inout) :: this
  type(c_ptr) :: ctxt
  integer(c_int) :: have_ctxt
  have_ctxt = atlas__connectivity__ctxt(this%CPTR_PGIBUG_A,ctxt)
  if( have_ctxt == 1 ) then
    call c_f_pointer(ctxt,this%access)
  else
    allocate( this%access )
    call atlas__connectivity__register_ctxt  ( this%CPTR_PGIBUG_A, c_loc(this%access) )
    call atlas__connectivity__register_update( this%CPTR_PGIBUG_A, c_funloc(update_access_c) )
    call atlas__connectivity__register_delete( this%CPTR_PGIBUG_A, c_funloc(delete_access_c) )

    this%access%connectivity_ptr = this%CPTR_PGIBUG_A
    call update_access(this%access)
  endif
end subroutine

subroutine update_access_c(this_ptr) bind(c)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  type(c_ptr), value :: this_ptr
  type(atlas_ConnectivityAccess), pointer :: this
  call c_f_pointer(this_ptr,this)
  call update_access(this)
end subroutine

subroutine update_access(this)
  use atlas_connectivity_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  type(atlas_ConnectivityAccess), intent(inout) :: this
  integer :: jrow

  type(c_ptr) :: values_cptr
  type(c_ptr) :: displs_cptr
  type(c_ptr) :: counts_cptr
  integer(ATLAS_KIND_IDX) :: values_size
  integer(ATLAS_KIND_IDX) :: displs_size
  integer(ATLAS_KIND_IDX) :: counts_size

  this%missing_value_ = atlas__Connectivity__missing_value(this%connectivity_ptr)
  call atlas__Connectivity__values(this%connectivity_ptr,values_cptr,values_size)
  call atlas__Connectivity__displs(this%connectivity_ptr,displs_cptr,displs_size)
  call atlas__Connectivity__counts(this%connectivity_ptr,counts_cptr,counts_size)

  call c_f_pointer(values_cptr, this%values_, (/values_size/))
  call c_f_pointer(displs_cptr, this%displs_, (/displs_size/))
  call c_f_pointer(counts_cptr, this%cols,    (/counts_size/))
  this%rows_   = atlas__Connectivity__rows(this%connectivity_ptr)
  if( associated( this%row ) ) deallocate(this%row)
  allocate( this%row(this%rows_) )
  this%maxcols_ = 0
  this%mincols_ = huge(this%maxcols_)
  do jrow=1,this%rows_
    this%row(jrow)%col => this%values_(this%displs_(jrow)+1:this%displs_(jrow+1)+1)
    this%maxcols_ = max(this%maxcols_, this%cols(jrow) )
    this%mincols_ = min(this%mincols_, this%cols(jrow) )
  enddo
  if( associated( this%padded_) ) call update_padded(this)
end subroutine

subroutine update_padded(this)
  class(atlas_ConnectivityAccess), intent(inout) :: this
  integer(ATLAS_KIND_IDX) :: jrow, jcol
  if( associated(this%padded_) ) call atlas_deallocate_managedmem( this%padded_ )
  call atlas_allocate_managedmem( this%padded_, [this%maxcols_,this%rows()] )
  if( associated(this%padded_) ) then
    this%padded_(:,:) = this%missing_value_
    do jrow=1,this%rows()
      do jcol=1,this%cols(jrow)
        this%padded_(jcol,jrow) = this%value(jcol,jrow)
      enddo
    enddo
  endif
end subroutine

subroutine delete_access_c(this_ptr) bind(c)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  type(c_ptr), value :: this_ptr
  type(atlas_ConnectivityAccess), pointer :: this
  call c_f_pointer(this_ptr,this)
  call delete_access(this)
end subroutine

subroutine delete_access(this)
  use, intrinsic :: iso_c_binding, only : c_int
  use atlas_connectivity_c_binding
  type(atlas_ConnectivityAccess), pointer, intent(inout) :: this
  if( associated( this%row ) )    deallocate(this%row)
  if( associated( this%padded_) ) call atlas_deallocate_managedmem(this%padded_)
#ifndef _CRAYFTN
  ! Cray compiler bug (CCE/8.7) leads to following runtime error here.
  !     lib-4412 : UNRECOVERABLE library error
  !       An argument in the DEALLOCATE statement is a disassociated pointer, an
  !       unallocated array, or a pointer not allocated as a pointer.
  ! --> Does this lead to a (tiny) memory leak?
  deallocate(this)
#endif
end subroutine

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Connectivity__final_auto(this)
  type(atlas_Connectivity), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Connectivity__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_MultiBlockConnectivity__final_auto(this)
  type(atlas_MultiBlockConnectivity), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_MultiBlockConnectivity__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!-------------------------------------------------------------------------------

end module atlas_connectivity_module

