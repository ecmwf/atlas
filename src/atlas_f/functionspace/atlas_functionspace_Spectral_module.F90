! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_functionspace_Spectral_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int
use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
#if ATLAS_HAVE_TRANS
use atlas_Trans_module, only: atlas_Trans
#endif
use atlas_Config_module, only: atlas_Config

implicit none

private :: c_ptr, c_int
private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Config
#if ATLAS_HAVE_TRANS
private :: atlas_Trans
#endif

public :: atlas_functionspace_Spectral

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_Spectral

! Purpose :
! -------
!   *atlas_functionspace_Spectral* : Interpretes spectral fields

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: gather_field
  procedure, private :: scatter_field
  procedure, private :: gather_fieldset
  procedure, private :: scatter_fieldset

  generic, public :: gather  =>  gather_field,  gather_fieldset
  generic, public :: scatter => scatter_field, scatter_fieldset

  procedure, private :: norm_scalar
  procedure, private :: norm_array
  generic, public :: norm => norm_scalar, norm_array

  procedure, public :: truncation
  procedure, public :: nb_spectral_coefficients
  procedure, public :: nb_spectral_coefficients_global
  procedure, public :: levels

  procedure, public :: nump
  procedure, public :: nmyms
  procedure, public :: nasm0
  procedure, public :: nvalue

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_functionspace_Spectral__final_auto
#endif

END TYPE atlas_functionspace_Spectral

interface atlas_functionspace_Spectral
  module procedure atlas_functionspace_Spectral__cptr
  module procedure atlas_functionspace_Spectral__config
#if ATLAS_HAVE_TRANS
  module procedure atlas_functionspace_Spectral__trans
#endif
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

function atlas_functionspace_Spectral__cptr(cptr) result(this)
  type(atlas_functionspace_Spectral) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_functionspace_Spectral__config(truncation,levels) result(this)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: this
  integer(c_int), intent(in)           :: truncation
  integer(c_int), intent(in), optional :: levels

  type(atlas_Config) :: options
  options = atlas_Config()

  call options%set("truncation",truncation)
  if( present(levels) ) call options%set("levels",levels)

  call this%reset_c_ptr( atlas__SpectralFunctionSpace__new__config(options%CPTR_PGIBUG_B) )
  call options%final()

  call this%return()
end function

#if ATLAS_HAVE_TRANS
function atlas_functionspace_Spectral__trans(trans,levels) result(this)
  use atlas_functionspace_spectral_c_binding
  type(atlas_functionspace_Spectral) :: this
  type(atlas_Trans), intent(in) :: trans
  integer(c_int), intent(in), optional :: levels

  type(atlas_Config) :: options
  options = atlas_Config()

  if( present(levels) ) call options%set("levels",levels)

  call this%reset_c_ptr( atlas__SpectralFunctionSpace__new__trans(trans%CPTR_PGIBUG_A, &
      options%CPTR_PGIBUG_B ) )
  call options%final()

  call this%return()
end function
#endif

subroutine gather_field(this,local,global)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__SpectralFunctionSpace__gather(this%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A)
end subroutine

subroutine scatter_field(this,global,local)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__SpectralFunctionSpace__scatter(this%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A)
end subroutine

subroutine gather_fieldset(this,local,global)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__SpectralFunctionSpace__gather_fieldset(this%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A)
end subroutine

subroutine scatter_fieldset(this,global,local)
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__SpectralFunctionSpace__scatter_fieldset(this%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A)
end subroutine

subroutine norm_scalar(this,field,norm,rank)
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  real(c_double), intent(out) :: norm
  integer(c_int), optional :: rank
  integer :: opt_rank
  real(c_double) :: norm_array(1)
  opt_rank = 0
  if( present(rank) ) opt_rank = rank
  call atlas__SpectralFunctionSpace__norm(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A,norm_array,opt_rank)
  norm = norm_array(1)
end subroutine

subroutine norm_array(this,field,norm,rank)
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use atlas_functionspace_spectral_c_binding
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  real(c_double), intent(inout) :: norm(:)
  integer(c_int), optional :: rank
  integer :: opt_rank
  opt_rank = 0
  if( present(rank) ) opt_rank = rank
  call atlas__SpectralFunctionSpace__norm(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A,norm,opt_rank)
end subroutine

function levels( this )
  use atlas_functionspace_spectral_c_binding
  integer :: levels
  class(atlas_functionspace_Spectral) :: this
  call atlas__SpectralFunctionSpace__levels( this%c_ptr(), levels )
end function

!-------------------------------------------------------------------------------
! Experimental

function truncation( this )
  use atlas_functionspace_spectral_c_binding
  integer :: truncation
  class(atlas_functionspace_Spectral) :: this
  call atlas__SpectralFunctionSpace__truncation( this%c_ptr(), truncation )
end function

function nb_spectral_coefficients( this )
  use atlas_functionspace_spectral_c_binding
  integer :: nb_spectral_coefficients
  class(atlas_functionspace_Spectral) :: this
  call atlas__SpectralFunctionSpace__nspec2( this%c_ptr(), nb_spectral_coefficients )
end function


function nb_spectral_coefficients_global( this )
  use atlas_functionspace_spectral_c_binding
  integer :: nb_spectral_coefficients_global
  class(atlas_functionspace_Spectral) :: this
  call atlas__SpectralFunctionSpace__nspec2g( this%c_ptr(), nb_spectral_coefficients_global )
end function


function nump( this )
  use atlas_functionspace_spectral_c_binding
  integer :: nump
  class(atlas_functionspace_Spectral) :: this
  call atlas__SpectralFunctionSpace__nump( this%c_ptr(), nump )
end function

function nmyms(this)
  use atlas_functionspace_spectral_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nmyms(:)
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(c_ptr) :: nmyms_c_ptr
  integer(c_int) :: size
  call atlas__SpectralFunctionSpace__nmyms(this%c_ptr(), nmyms_c_ptr, size)
  call c_f_pointer ( nmyms_c_ptr , nmyms , (/size/) )
end function

function nasm0(this)
  use atlas_functionspace_spectral_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nasm0(:)
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(c_ptr) :: nasm0_c_ptr
  integer(c_int), pointer :: nasm0_f_ptr(:)
  integer(c_int) :: size
  call atlas__SpectralFunctionSpace__nasm0(this%c_ptr(), nasm0_c_ptr, size)
  call c_f_pointer ( nasm0_c_ptr , nasm0_f_ptr , (/size/) )
  nasm0(0:) => nasm0_f_ptr(:)
end function

function nvalue(this)
  use atlas_functionspace_spectral_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nvalue(:)
  class(atlas_functionspace_Spectral), intent(in) :: this
  type(c_ptr) :: nvalue_c_ptr
  integer(c_int) :: size
  call atlas__SpectralFunctionSpace__nvalue(this%c_ptr(), nvalue_c_ptr, size)
  call c_f_pointer ( nvalue_c_ptr , nvalue , (/size/) )
end function




!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_functionspace_Spectral__final_auto(this)
  type(atlas_functionspace_Spectral), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_functionspace_Spectral__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_functionspace_Spectral_module

