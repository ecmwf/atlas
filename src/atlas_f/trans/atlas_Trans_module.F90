#include "atlas/atlas_f.h"

module atlas_Trans_module


use fckit_object_module, only: fckit_object
use fckit_owned_object_module, only: fckit_owned_object
use atlas_config_module, only : atlas_Config
use atlas_field_module, only : atlas_Field
use atlas_fieldset_module, only : atlas_FieldSet
use atlas_grid_module, only : atlas_Grid

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
  procedure :: handle
  procedure :: grid
  procedure :: truncation
  procedure :: nb_spectral_coefficients
  procedure :: nb_spectral_coefficients_global
  procedure :: nb_gridpoints
  procedure :: nb_gridpoints_global

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

  procedure, private :: gathspec_r1
  procedure, private :: gathspec_r2
  generic, public :: gathspec => gathspec_r1, gathspec_r2

  procedure, private :: specnorm_r1_scalar
  procedure, private :: specnorm_r2
  generic, public :: specnorm => specnorm_r1_scalar, specnorm_r2

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

#define THROW_ERROR call te("atlas_Trans_module.F90",__LINE__)

subroutine te(file,line)
  use fckit_exception_module, only : fckit_exception
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  call fckit_exception%throw( "Cannot use atlas_Trans since atlas is compiled without" // &
    & "ENABLE_TRANS=ON", file, line )
end subroutine

function atlas_Trans__ctor( grid, nsmax ) result(this)
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  use atlas_trans_c_binding
  type(atlas_Trans) :: this
  class(atlas_Grid), intent(in) :: grid
  integer, intent(in), optional :: nsmax
#if ATLAS_HAVE_TRANS
  if( present(nsmax) ) then
    call this%reset_c_ptr( atlas__Trans__new( grid%CPTR_PGIBUG_A, nsmax ) )
  else
    call this%reset_c_ptr( atlas__Trans__new( grid%CPTR_PGIBUG_A, 0 ) )
  endif
#else
  ! IGNORE
  call this%reset_c_ptr( c_null_ptr )
  FCKIT_SUPPRESS_UNUSED( grid )
  FCKIT_SUPPRESS_UNUSED( nsmax )
#endif
  call this%return()
end function atlas_Trans__ctor

function handle( this )
  use atlas_trans_c_binding
  integer :: handle
  class(atlas_Trans) :: this
#if ATLAS_HAVE_TRANS
  handle = atlas__Trans__handle (this%CPTR_PGIBUG_A)
#else
  THROW_ERROR
  handle = 0
  FCKIT_SUPPRESS_UNUSED( this )
#endif
end function

function truncation( this )
  use atlas_trans_c_binding
  integer :: truncation
  class(atlas_Trans) :: this
#if ATLAS_HAVE_TRANS
  truncation = atlas__Trans__truncation (this%CPTR_PGIBUG_A)
#else
  THROW_ERROR
  truncation = 0
  FCKIT_SUPPRESS_UNUSED( this )
#endif
end function

function nb_spectral_coefficients( this )
  use atlas_trans_c_binding
  integer :: nb_spectral_coefficients
  class(atlas_Trans) :: this
#if ATLAS_HAVE_TRANS
  nb_spectral_coefficients = atlas__Trans__nspec2 (this%CPTR_PGIBUG_A)
#else
  THROW_ERROR
  nb_spectral_coefficients = 0
  FCKIT_SUPPRESS_UNUSED( this )
#endif
end function

function nb_spectral_coefficients_global( this )
  use atlas_trans_c_binding
  integer :: nb_spectral_coefficients_global
  class(atlas_Trans) :: this
#if ATLAS_HAVE_TRANS
  nb_spectral_coefficients_global = atlas__Trans__nspec2g (this%CPTR_PGIBUG_A)
#else
  THROW_ERROR
  nb_spectral_coefficients_global = 0
  FCKIT_SUPPRESS_UNUSED( this )
#endif
end function

function nb_gridpoints( this )
  use atlas_trans_c_binding
  integer :: nb_gridpoints
  class(atlas_Trans) :: this
#if ATLAS_HAVE_TRANS
  nb_gridpoints = atlas__Trans__ngptot (this%CPTR_PGIBUG_A)
#else
  THROW_ERROR
  nb_gridpoints = 0
  FCKIT_SUPPRESS_UNUSED( this )
#endif
end function

function nb_gridpoints_global( this )
  use atlas_trans_c_binding
  integer :: nb_gridpoints_global
  class(atlas_Trans) :: this
#if ATLAS_HAVE_TRANS
  nb_gridpoints_global = atlas__Trans__ngptotg (this%CPTR_PGIBUG_A)
#else
  THROW_ERROR
  nb_gridpoints_global = 0
  FCKIT_SUPPRESS_UNUSED( this )
#endif
end function

function grid( this )
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  use atlas_trans_c_binding
  class(atlas_Trans) :: this
  type(atlas_Grid) :: grid
#if ATLAS_HAVE_TRANS
  grid = atlas_Grid( atlas__Trans__grid(this%CPTR_PGIBUG_A) )
  call grid%return()
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  grid = atlas_Grid( c_null_ptr )
  FCKIT_SUPPRESS_UNUSED( grid )
#endif
end function


subroutine dirtrans_fieldset(this, gpfields, spfields, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: gpfields
  class(atlas_FieldSet), intent(inout) :: spfields
  class(atlas_Config), intent(in), optional  :: config
#if ATLAS_HAVE_TRANS
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
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( gpfields )
  FCKIT_SUPPRESS_UNUSED( spfields )
  FCKIT_SUPPRESS_UNUSED( config )
#endif
end subroutine dirtrans_fieldset


subroutine invtrans_fieldset(this, spfields, gpfields, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: spfields
  class(atlas_FieldSet), intent(inout) :: gpfields
  class(atlas_Config), intent(in), optional  :: config
#if ATLAS_HAVE_TRANS
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
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( spfields )
  FCKIT_SUPPRESS_UNUSED( gpfields )
  FCKIT_SUPPRESS_UNUSED( config )
#endif
end subroutine invtrans_fieldset

subroutine dirtrans_field(this, gpfield, spfield, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: gpfield
  class(atlas_Field), intent(inout) :: spfield
  class(atlas_Config), intent(in), optional  :: config
#if ATLAS_HAVE_TRANS
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
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( gpfield )
  FCKIT_SUPPRESS_UNUSED( spfield )
  FCKIT_SUPPRESS_UNUSED( config )
#endif
end subroutine dirtrans_field

subroutine dirtrans_wind2vordiv_field(this, gpwind, spvor, spdiv, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  type(atlas_Field), intent(in)  :: gpwind
  type(atlas_Field), intent(inout) :: spvor
  type(atlas_Field), intent(inout) :: spdiv
  type(atlas_Config), intent(in), optional  :: config
#if ATLAS_HAVE_TRANS
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
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( gpwind )
  FCKIT_SUPPRESS_UNUSED( spvor )
  FCKIT_SUPPRESS_UNUSED( spdiv )
  FCKIT_SUPPRESS_UNUSED( config )
#endif

end subroutine dirtrans_wind2vordiv_field


subroutine invtrans_field(this, spfield, gpfield, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_Field), intent(inout) :: gpfield
  class(atlas_Config), intent(in), optional  :: config
#if ATLAS_HAVE_TRANS
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
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( spfield )
  FCKIT_SUPPRESS_UNUSED( gpfield )
  FCKIT_SUPPRESS_UNUSED( config )
#endif
end subroutine invtrans_field


subroutine invtrans_vordiv2wind_field(this, spvor, spdiv, gpwind, config)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spvor
  class(atlas_Field), intent(in)  :: spdiv
  class(atlas_Field), intent(inout) :: gpwind
  class(atlas_Config), intent(in), optional  :: config
#if ATLAS_HAVE_TRANS
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
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( spvor )
  FCKIT_SUPPRESS_UNUSED( spdiv )
  FCKIT_SUPPRESS_UNUSED( gpwind )
  FCKIT_SUPPRESS_UNUSED( config )
#endif

end subroutine invtrans_vordiv2wind_field


subroutine invtrans_grad_field(this, spfield, gpfield)
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_Field), intent(inout) :: gpfield
#if ATLAS_HAVE_TRANS
  type(atlas_Config) :: config
  config = atlas_Config()
  call atlas__Trans__invtrans_grad_field( this%CPTR_PGIBUG_A, &
    &                          spfield%CPTR_PGIBUG_A, &
    &                          gpfield%CPTR_PGIBUG_A, &
    &                          config%CPTR_PGIBUG_B)
  call config%final()
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( spfield )
  FCKIT_SUPPRESS_UNUSED( gpfield )
#endif
end subroutine invtrans_grad_field



subroutine gathspec_r1(this, local, global)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: local(:)
  real(c_double), intent(inout) :: global(:)
#if ATLAS_HAVE_TRANS
  call atlas__Trans__gathspec(this%CPTR_PGIBUG_A, 1, (/1/), local, global )
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( local )
  FCKIT_SUPPRESS_UNUSED( global )
#endif
end subroutine gathspec_r1

subroutine gathspec_r2(this, local, global)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_array_module, only: array_view1d
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: local(:,:)
  real(c_double), intent(inout) :: global(:,:)
#if ATLAS_HAVE_TRANS
  real(c_double), pointer :: local_view(:), global_view(:)
  integer :: destination(size(local,1))
  destination(:) = 1
  local_view => array_view1d(local)
  global_view => array_view1d(global)
  call atlas__Trans__gathspec(this%CPTR_PGIBUG_A, size(local,1), destination, local_view, global_view )
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( local )
  FCKIT_SUPPRESS_UNUSED( global )
#endif
end subroutine gathspec_r2


subroutine specnorm_r1_scalar(this, spectra, norm, rank)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: spectra(:)
  real(c_double), intent(out) :: norm
  integer, optional :: rank ! MPI rank
#if ATLAS_HAVE_TRANS
  integer :: rank_opt
  real(c_double) :: norms(1)
  rank_opt = 0
  if( present(rank) ) rank_opt = rank
  call atlas__Trans__specnorm(this%CPTR_PGIBUG_A, 1, spectra, norms, rank_opt )
  norm = norms(1)
#else
  norm=0
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( spectra )
  FCKIT_SUPPRESS_UNUSED( rank )
#endif
end subroutine

subroutine specnorm_r2(this, spectra, norm, rank)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_array_module, only: array_view1d
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: spectra(:,:)
  real(c_double), intent(inout) :: norm(:)
  integer, optional :: rank ! MPI rank
#if ATLAS_HAVE_TRANS
  integer :: rank_opt
  real(c_double), pointer :: spectra_view(:)
  rank_opt = 0
  if( present(rank) ) rank_opt = rank
  spectra_view => array_view1d(spectra)
  call atlas__Trans__specnorm(this%CPTR_PGIBUG_A, size(spectra,1), spectra_view, norm, rank_opt )
#else
  THROW_ERROR
  FCKIT_SUPPRESS_UNUSED( this )
  FCKIT_SUPPRESS_UNUSED( spectra )
  FCKIT_SUPPRESS_UNUSED( norm )
  FCKIT_SUPPRESS_UNUSED( rank )
#endif
end subroutine

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
