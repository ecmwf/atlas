! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Trans_module


use fckit_object_module, only: fckit_object
use fckit_owned_object_module, only: fckit_owned_object
use atlas_config_module, only : atlas_Config
use atlas_field_module, only : atlas_Field
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_grid_module, only : atlas_Grid
use atlas_functionspace_module, only : atlas_functionspace

implicit none

public :: atlas_Trans

private

!-----------------------------
! atlas_Trans                 !
!-----------------------------

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Trans

! Purpose :
! -------
!   *Trans* : To do spectral transforms

! Methods :
! -------

! Author :
! ------
!   12-Mar-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, nopass :: has_backend
  procedure, nopass :: set_backend
  procedure, nopass :: backend

  procedure :: handle
  procedure :: grid
  procedure :: truncation
  procedure :: spectral
#if 0
  procedure :: nb_spectral_coefficients
  procedure :: nb_spectral_coefficients_global
  procedure :: nb_gridpoints
  procedure :: nb_gridpoints_global
#endif

  procedure, private :: dirtrans_field
  procedure, private :: dirtrans_fieldset
  procedure, public :: dirtrans_wind2vordiv => dirtrans_wind2vordiv_field
  generic, public :: dirtrans => &
    & dirtrans_field, &
    & dirtrans_fieldset

  procedure, private :: invtrans_field
  procedure, private :: invtrans_fieldset
  procedure, public :: invtrans_vordiv2wind => invtrans_vordiv2wind_field
  generic, public :: invtrans => &
    & invtrans_field, &
    & invtrans_fieldset

  procedure, private :: invtrans_grad_field
  generic, public :: invtrans_grad => &
    & invtrans_grad_field

! removed
#if 0
  procedure, private :: gathspec_r1
  procedure, private :: gathspec_r2
  generic, public :: gathspec => gathspec_r1, gathspec_r2

  procedure, private :: specnorm_r1_scalar
  procedure, private :: specnorm_r2
  generic, public :: specnorm => specnorm_r1_scalar, specnorm_r2
#endif

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Trans__final_auto
#endif

END TYPE atlas_Trans

!------------------------------------------------------------------------------

interface atlas_Trans
  module procedure atlas_Trans__ctor
end interface

!------------------------------------------------------------------------------

private :: fckit_owned_object
private :: fckit_object
private :: atlas_Config
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Grid

!========================================================
contains
!========================================================

! -----------------------------------------------------------------------------
! Trans routines

function backend()
  use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
  use atlas_trans_c_binding
  use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
  character(len=:), allocatable :: backend
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: value_size
  call atlas__Trans__backend(value_cptr,value_size)
  allocate( character(len=value_size) :: backend )
  backend = c_ptr_to_string(value_cptr)
  call c_ptr_free(value_cptr)
end function

function has_backend( backend )
  use, intrinsic :: iso_c_binding, only: c_int
  use fckit_c_interop_module, only : c_str
  use atlas_trans_c_binding
  logical :: has_backend
  character(len=*) :: backend
  integer(c_int) :: has_backend_int
  has_backend_int =  atlas__Trans__has_backend( c_str(backend) )
  if( has_backend_int == 1 ) then
    has_backend = .True.
  else
    has_backend = .False.
  end if
end function

subroutine set_backend( backend )
  use atlas_trans_c_binding
  use fckit_c_interop_module, only : c_str
  character(len=*), intent(in) :: backend
  call atlas__Trans__set_backend( c_str(backend) )
end subroutine


function atlas_Trans__ctor( grid, nsmax, config ) result(this)
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  use atlas_trans_c_binding
  type(atlas_Trans) :: this
  class(atlas_Grid), intent(in) :: grid
  integer, intent(in), optional :: nsmax
  type(atlas_Config), intent(in), optional :: config
  if( present( config ) ) then
    if( present(nsmax) ) then
      call this%reset_c_ptr( atlas__Trans__new_config( grid%CPTR_PGIBUG_A, nsmax, config%CPTR_PGIBUG_B ) )
    else
      call this%reset_c_ptr( atlas__Trans__new_config( grid%CPTR_PGIBUG_A, 0, config%CPTR_PGIBUG_B ) )
    endif
  else
    if( present(nsmax) ) then
      call this%reset_c_ptr( atlas__Trans__new( grid%CPTR_PGIBUG_A, nsmax ) )
    else
      call this%reset_c_ptr( atlas__Trans__new( grid%CPTR_PGIBUG_A, 0 ) )
    endif
  endif
  call this%return()
end function atlas_Trans__ctor

function handle( this )
  use atlas_trans_c_binding
  integer :: handle
  class(atlas_Trans) :: this
  handle = atlas__Trans__handle (this%CPTR_PGIBUG_A)
end function

function truncation( this )
  use atlas_trans_c_binding
  integer :: truncation
  class(atlas_Trans) :: this
  truncation = atlas__Trans__truncation (this%CPTR_PGIBUG_A)
end function

function spectral( this )
  use atlas_trans_c_binding
  type(atlas_FunctionSpace) :: spectral
  class(atlas_Trans) :: this
  spectral = atlas_FunctionSpace( atlas__Trans__spectral (this%CPTR_PGIBUG_A) )
  call spectral%return()
end function

function grid( this )
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  use atlas_trans_c_binding
  class(atlas_Trans) :: this
  type(atlas_Grid) :: grid
  grid = atlas_Grid( atlas__Trans__grid(this%CPTR_PGIBUG_A) )
  call grid%return()
end function


subroutine dirtrans_fieldset(this, gpfields, spfields, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: gpfields
  class(atlas_FieldSet), intent(inout) :: spfields
  class(atlas_Config), intent(in), optional  :: config
  type(atlas_Config) :: p

  if( present(config) ) then
    call p%reset_c_ptr( config%CPTR_PGIBUG_B )
  else
    p = atlas_Config()
  endif

  call atlas__Trans__dirtrans_fieldset( this%CPTR_PGIBUG_A,     &
    &                          gpfields%CPTR_PGIBUG_A, &
    &                          spfields%CPTR_PGIBUG_A, &
    &                          p%CPTR_PGIBUG_B )

  if( .not. present(config) ) then
    call p%final()
  endif
end subroutine dirtrans_fieldset


subroutine invtrans_fieldset(this, spfields, gpfields, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: spfields
  class(atlas_FieldSet), intent(inout) :: gpfields
  class(atlas_Config), intent(in), optional  :: config
  type(atlas_Config) :: p

  if( present(config) ) then
    call p%reset_c_ptr( config%CPTR_PGIBUG_B )
  else
    p = atlas_Config()
  endif

  call atlas__Trans__invtrans_fieldset( this%CPTR_PGIBUG_A,     &
    &                          spfields%CPTR_PGIBUG_A, &
    &                          gpfields%CPTR_PGIBUG_A, &
    &                          p%CPTR_PGIBUG_B )

  if( .not. present(config) ) then
    call p%final()
  endif
end subroutine invtrans_fieldset

subroutine dirtrans_field(this, gpfield, spfield, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: gpfield
  class(atlas_Field), intent(inout) :: spfield
  class(atlas_Config), intent(in), optional  :: config
  type(atlas_Config) :: p

  if( present(config) ) then
    call p%reset_c_ptr( config%CPTR_PGIBUG_B )
  else
    p = atlas_Config()
  endif

  call atlas__Trans__dirtrans_field( this%CPTR_PGIBUG_A, &
    &                          gpfield%CPTR_PGIBUG_A, &
    &                          spfield%CPTR_PGIBUG_A, &
    &                          p%CPTR_PGIBUG_B )

  if( .not. present(config) ) then
    call p%final()
  endif
end subroutine dirtrans_field

subroutine dirtrans_wind2vordiv_field(this, gpwind, spvor, spdiv, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  type(atlas_Field), intent(in)  :: gpwind
  type(atlas_Field), intent(inout) :: spvor
  type(atlas_Field), intent(inout) :: spdiv
  type(atlas_Config), intent(in), optional  :: config
  type(atlas_Config) :: p

  if( present(config) ) then
    call p%reset_c_ptr( config%CPTR_PGIBUG_B )
  else
    p = atlas_Config()
  endif

  call atlas__Trans__dirtrans_wind2vordiv_field( this%CPTR_PGIBUG_A, &
    &                          gpwind%CPTR_PGIBUG_A, &
    &                          spvor%CPTR_PGIBUG_A, &
    &                          spdiv%CPTR_PGIBUG_A, &
    &                          p%CPTR_PGIBUG_B )

  if( .not. present(config) ) then
    call p%final()
  endif

end subroutine dirtrans_wind2vordiv_field


subroutine invtrans_field(this, spfield, gpfield, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_Field), intent(inout) :: gpfield
  class(atlas_Config), intent(in), optional  :: config
  type(atlas_Config) :: p

  if( present(config) ) then
    call p%reset_c_ptr( config%CPTR_PGIBUG_B )
  else
    p = atlas_Config()
  endif

  call atlas__Trans__invtrans_field( this%CPTR_PGIBUG_A, &
    &                          spfield%CPTR_PGIBUG_A, &
    &                          gpfield%CPTR_PGIBUG_A, &
    &                          p%CPTR_PGIBUG_B )

  if( .not. present(config) ) then
    call p%final()
  endif
end subroutine invtrans_field


subroutine invtrans_vordiv2wind_field(this, spvor, spdiv, gpwind, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spvor
  class(atlas_Field), intent(in)  :: spdiv
  class(atlas_Field), intent(inout) :: gpwind
  class(atlas_Config), intent(in), optional  :: config
  type(atlas_Config) :: p

  if( present(config) ) then
    call p%reset_c_ptr( config%CPTR_PGIBUG_B )
  else
    p = atlas_Config()
  endif

  call atlas__Trans__invtrans_vordiv2wind_field( this%CPTR_PGIBUG_A, &
    &                          spvor%CPTR_PGIBUG_A, &
    &                          spdiv%CPTR_PGIBUG_A, &
    &                          gpwind%CPTR_PGIBUG_A, &
    &                          p%CPTR_PGIBUG_B )

  if( .not. present(config) ) then
    call p%final()
  endif

end subroutine invtrans_vordiv2wind_field


subroutine invtrans_grad_field(this, spfield, gpfield)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_Field), intent(inout) :: gpfield
  type(atlas_Config) :: config
  config = atlas_Config()
  call atlas__Trans__invtrans_grad_field( this%CPTR_PGIBUG_A, &
    &                          spfield%CPTR_PGIBUG_A, &
    &                          gpfield%CPTR_PGIBUG_A, &
    &                          config%CPTR_PGIBUG_B)
  call config%final()
end subroutine invtrans_grad_field

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_Trans__final_auto(this)
  type(atlas_Trans), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Trans__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_Trans_module
