! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Method_module

use fckit_owned_object_module, only : fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_Method

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Method

! Purpose :
! -------
!   *Method* :
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done

! Methods :
! -------
!   name : The name or tag this function space was created with

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: name => atlas_Method__name
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Method__final_auto
#endif
END TYPE atlas_Method

interface atlas_Method
  module procedure atlas_Method__cptr
end interface

!========================================================
contains
!========================================================

function atlas_Method__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Method) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_Method__name(this) result(name)
  use atlas_Method_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Method), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__Method__name(this%CPTR_PGIBUG_A)
  name = c_ptr_to_string(name_c_str)
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Method__final_auto(this)
  type(atlas_Method), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Method__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_Method_module

