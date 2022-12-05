! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Trace_module

use fckit_shared_object_module, only : fckit_shared_object

implicit none

private :: fckit_shared_object

public :: atlas_Trace


private

!-----------------------------
! atlas_Trace                !
!-----------------------------

type, extends(fckit_shared_object) :: atlas_Trace
contains
! Public methods

  procedure, public :: running
  procedure, public :: start
  procedure, public :: stop
  procedure, public :: pause
  procedure, public :: resume
  procedure, public :: elapsed

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Trace__final_auto
#endif
end type

interface atlas_Trace
  module procedure atlas_Trace__loc
  module procedure atlas_Trace__labels_1
  module procedure atlas_Trace__labels_2
  module procedure atlas_Trace__labels_3
  module procedure atlas_Trace__labels_4
  module procedure atlas_Trace__labels_5
end interface

!========================================================
contains
!========================================================

function atlas_Trace__loc(file,line,title) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  use fckit_c_interop_module
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  call this%reset_c_ptr( new_atlas_Trace( c_str(file), line, c_str(title) ), fckit_c_deleter(delete_atlas_Trace) )
  call this%return()
end function

function atlas_Trace__labels_1(file,line,title,label) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  use fckit_c_interop_module
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label
  call this%reset_c_ptr( new_atlas_Trace_labels_1( c_str(file), line, c_str(title), c_str(label) ), &
    &  fckit_c_deleter(delete_atlas_Trace) )
  call this%return()
end function
function atlas_Trace__labels_2(file,line,title,label1,label2) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  use fckit_c_interop_module
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  call this%reset_c_ptr( new_atlas_Trace_labels_2( c_str(file), line, c_str(title), c_str(label1), c_str(label2) ), &
    &  fckit_c_deleter(delete_atlas_Trace) )
  call this%return()
end function
function atlas_Trace__labels_3(file,line,title,label1,label2,label3) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  use fckit_c_interop_module
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  character(len=*) , intent(in) :: label3
  call this%reset_c_ptr( new_atlas_Trace_labels_3( c_str(file), line, c_str(title), c_str(label1), c_str(label2), &
    & c_str(label3) ), fckit_c_deleter(delete_atlas_Trace) )
  call this%return()
end function
function atlas_Trace__labels_4(file,line,title,label1,label2,label3,label4) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  use fckit_c_interop_module
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  character(len=*) , intent(in) :: label3
  character(len=*) , intent(in) :: label4
  call this%reset_c_ptr( new_atlas_Trace_labels_4( c_str(file), line, c_str(title), c_str(label1), c_str(label2), &
    & c_str(label3), c_str(label4) ), fckit_c_deleter(delete_atlas_Trace) )
  call this%return()
end function
function atlas_Trace__labels_5(file,line,title,label1,label2,label3,label4,label5) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Trace_c_binding
  use fckit_c_interop_module
  type(atlas_Trace) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label1
  character(len=*) , intent(in) :: label2
  character(len=*) , intent(in) :: label3
  character(len=*) , intent(in) :: label4
  character(len=*) , intent(in) :: label5
  call this%reset_c_ptr( new_atlas_Trace_labels_5( c_str(file), line, c_str(title), c_str(label1), c_str(label2), &
    & c_str(label3), c_str(label4), c_str(label5) ), fckit_c_deleter(delete_atlas_Trace) )
  call this%return()
end function
!-------------------------------------------------------------------------------

function running( this )
  use atlas_Trace_c_binding
  logical :: running
  class(atlas_Trace) :: this
  if( atlas_Trace__running( this%CPTR_PGIBUG_B ) == 0 ) then
    running = .False.
  else
    running = .True.
  endif
end function

!-------------------------------------------------------------------------------

function elapsed( this )
  use, intrinsic :: iso_c_binding, only : c_double
  use atlas_Trace_c_binding
  real(c_double) :: elapsed
  class(atlas_Trace) :: this
  elapsed = atlas_Trace__elapsed( this%CPTR_PGIBUG_B )
end function

!-------------------------------------------------------------------------------

subroutine start( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__start( this%CPTR_PGIBUG_B )
end subroutine

!-------------------------------------------------------------------------------

subroutine stop( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__stop( this%CPTR_PGIBUG_B )
end subroutine

!-------------------------------------------------------------------------------

subroutine pause( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__pause( this%CPTR_PGIBUG_B )
end subroutine

!-------------------------------------------------------------------------------

subroutine resume( this )
  use atlas_Trace_c_binding
  class(atlas_Trace) :: this
  call atlas_Trace__resume( this%CPTR_PGIBUG_B )
end subroutine

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Trace__final_auto(this)
  type(atlas_Trace), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Trace__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_Trace_module

