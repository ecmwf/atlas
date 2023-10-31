! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Redistribution_module

use, intrinsic :: iso_c_binding, only : c_ptr
use atlas_config_module, only : atlas_Config
use atlas_functionspace_module, only : atlas_FunctionSpace
use fckit_owned_object_module, only: fckit_owned_object

implicit none

public :: atlas_Redistribution

private

!-----------------------------
! atlas_Redistribution     !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Redistribution

! Purpose :
! -------
!   *Redistribution* : Object passed from atlas to inspect redistribution

! Methods :
! -------

! Author :
! ------

!------------------------------------------------------------------------------
contains

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Redistribution__final_auto
#endif
END TYPE atlas_Redistribution

!------------------------------------------------------------------------------

interface atlas_Redistribution
  module procedure ctor_cptr
  module procedure ctor_create
end interface

private :: c_ptr
private :: fckit_owned_object

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! Redistribution routines

function empty_config() result(config)
  type(atlas_Config) :: config
  config = atlas_Config()
  call config%return()
end function

function ctor_cptr( cptr ) result(this)
  use atlas_redistribution_c_binding
  type(atlas_Redistribution) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function ctor_create( fspace1, fspace2, redist_name ) result(this)
  use atlas_redistribution_c_binding
  class(atlas_FunctionSpace), intent(in) :: fspace1, fspace2
  character(len=*), intent(in), optional :: redist_name
  type(atlas_Redistribution) :: this
  type(atlas_Config) :: config
  config = empty_config()
  call config%set("type", "RedistributeGeneric")
  if (present(redist_name)) call config%set("type", redist_name)
  call this%reset_c_ptr( atlas__Redistribution__new__config(fspace1%CPTR_PGIBUG_A, fspace2%CPTR_PGIBUG_A, config%CPTR_PGIBUG_B) )
  call config%final()
  call this%return()
end function

! ----------------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Redistribution__final_auto(this)
  type(atlas_Redistribution), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Redistribution__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! ----------------------------------------------------------------------------------------

end module atlas_Redistribution_module
