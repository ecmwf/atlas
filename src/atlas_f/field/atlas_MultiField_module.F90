! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_multifield_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_Config_module, only: atlas_Config
use atlas_fieldset_module, only: atlas_fieldset

implicit none

public :: atlas_MultiField

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_MultiField

! Purpose :
! -------
!   *MultiField* : Object containing field data and Metadata

! Methods :
! -------
!   name : The name or tag this field was created with
!   data : Return the field as a fortran array of specified shape
!   Metadata : Return object that can contain a variety of Metadata

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*
!   29-Aug-2023 Slavko Brdar         *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public  :: MultiField__fieldset
  procedure, public  :: MultiField__size
  !procedure, public  :: MultiField__data
  generic :: fieldset => MultiField__fieldset
  generic :: size => MultiField__size

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_MultiField__final_auto
#endif

END TYPE

interface atlas_MultiField
  module procedure atlas_MultiField__cptr
  module procedure atlas_MultiField__create
end interface

private :: fckit_owned_object
private :: atlas_Config

!========================================================
contains
!========================================================

!-------------------------------------------------------------------------------

function atlas_MultiField__cptr(cptr) result(field)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_MultiField) :: field
  type(c_ptr), intent(in) :: cptr
  call field%reset_c_ptr( cptr )
  call field%return()
end function

!-------------------------------------------------------------------------------

function atlas_MultiField__create(params) result(field)
  use atlas_multifield_c_binding
  type(atlas_MultiField) :: field
  class(atlas_Config), intent(in) :: params
  field = atlas_MultiField__cptr( atlas__MultiField__create(params%CPTR_PGIBUG_B) )
  call field%return()
end function

!-------------------------------------------------------------------------------

function MultiField__size(this) result(size)
  use atlas_multifield_c_binding
  class(atlas_MultiField), intent(in) :: this
  integer :: size
  size = atlas__MultiField__size(this%CPTR_PGIBUG_B)
end function

!-------------------------------------------------------------------------------
function MultiField__fieldset(this) result(fset)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_multifield_c_binding
  class(atlas_MultiField), intent(in) :: this
  type(c_ptr) :: fset_cptr
  type(atlas_FieldSet) :: fset
  fset_cptr = atlas__MultiField__fieldset(this%CPTR_PGIBUG_B)
  fset = atlas_FieldSet( fset_cptr )
  call fset%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_MultiField__final_auto(this)
  type(atlas_MultiField), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_MultiField__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!-------------------------------------------------------------------------------

end module atlas_multifield_module

