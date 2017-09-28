
#include "atlas/atlas_f.h"

module atlas_Trans_module


use fckit_object_module, only: fckit_object
use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted
private :: fckit_object

public :: atlas_Trans
public :: atlas_TransParameters

private

!-----------------------------
! atlas_Trans                 !
!-----------------------------

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_Trans

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

  procedure, private :: dirtrans_field_nodes
  procedure, private :: dirtrans_field
  procedure, private :: dirtrans_fieldset_nodes
  procedure, private :: dirtrans_fieldset
  procedure, public :: dirtrans_wind2vordiv => dirtrans_wind2vordiv_field
  generic, public :: dirtrans => &
    & dirtrans_field, &
    & dirtrans_fieldset, &
    & dirtrans_fieldset_nodes, &
    & dirtrans_field_nodes

  procedure, private :: invtrans_field_nodes
  procedure, private :: invtrans_field
  procedure, private :: invtrans_fieldset_nodes
  procedure, private :: invtrans_fieldset
  procedure, public :: invtrans_vordiv2wind => invtrans_vordiv2wind_field
  generic, public :: invtrans => &
    & invtrans_field, &
    & invtrans_fieldset, &
    & invtrans_field_nodes, &
    & invtrans_fieldset_nodes

  procedure, private :: invtrans_grad_field_nodes
  generic, public :: invtrans_grad => &
    & invtrans_grad_field_nodes

  procedure, private :: gathspec_r1
  procedure, private :: gathspec_r2
  generic, public :: gathspec => gathspec_r1, gathspec_r2

  procedure, private :: specnorm_r1_scalar
  procedure, private :: specnorm_r2
  generic, public :: specnorm => specnorm_r1_scalar, specnorm_r2

  procedure, public :: delete => atlas_Trans__delete
  procedure, public :: copy => atlas_Trans__copy

END TYPE atlas_Trans

!------------------------------------------------------------------------------

interface atlas_Trans
  module procedure atlas_Trans__ctor
end interface

!------------------------------------------------------------------------------

TYPE, extends(fckit_object) :: atlas_TransParameters

! Purpose :
! -------
!   *TransParameters* : Extra information to pass to dirtrans and invtrans

! Author :
! ------
!   20-Mar-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, public :: delete => atlas_TransParameters__delete
  procedure, public :: copy => atlas_TransParameters__copy

END TYPE atlas_TransParameters

!------------------------------------------------------------------------------

interface atlas_TransParameters
  module procedure atlas_TransParameters__ctor
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! Trans routines

#define THROW_ERROR call te("atlas_Trans_module.F90",__LINE__)

subroutine te(file,line)
  use atlas_Error_module, only: atlas_code_location, atlas_throw_usererror
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  call atlas_throw_usererror("Cannot use atlas_Trans since atlas is compiled without" // &
    & "ENABLE_TRANS=ON",atlas_code_location(file,line))
end subroutine

function atlas_Trans__ctor( grid, nsmax ) result(trans)
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  use atlas_trans_c_binding
  use atlas_Grid_module, only: atlas_Grid
  type(atlas_Trans) :: trans
  class(atlas_Grid), intent(in) :: grid
  integer, intent(in), optional :: nsmax
#ifdef ATLAS_HAVE_TRANS
  if( present(nsmax) ) then
    call trans%reset_c_ptr( atlas__Trans__new( grid%c_ptr(), nsmax ) )
  else
    call trans%reset_c_ptr( atlas__Trans__new( grid%c_ptr(), 0 ) )
  endif
#else
  ! IGNORE
  call trans%reset_c_ptr( c_null_ptr )
#endif
end function atlas_Trans__ctor

function atlas_TransParameters__ctor() result(params)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  type(atlas_TransParameters) :: params
#ifdef ATLAS_HAVE_TRANS
  call params%reset_c_ptr( atlas__TransParameters__new() )
#else
  ! IGNORE
  call params%reset_c_ptr( c_null_ptr)
#endif
end function atlas_TransParameters__ctor

subroutine atlas_Trans__delete( this )
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  class(atlas_Trans), intent(inout) :: this
#ifdef ATLAS_HAVE_TRANS
  call atlas__Trans__delete(this%c_ptr());
#else
  ! IGNORE
  call this%reset_c_ptr( c_null_ptr )
#endif
end subroutine


subroutine atlas_Trans__copy(this,obj_in)
  class(atlas_Trans), intent(inout) :: this
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine



subroutine atlas_TransParameters__delete( this )
  use atlas_trans_c_binding
  class(atlas_TransParameters), intent(inout) :: this
#ifdef ATLAS_HAVE_TRANS
  call atlas__TransParameters__delete(this%c_ptr());
#else
  ! IGNORE
#endif
end subroutine


subroutine atlas_TransParameters__copy(this,obj_in)
  class(atlas_TransParameters), intent(inout) :: this
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine

function handle( this )
  use atlas_trans_c_binding
  integer :: handle
  class(atlas_Trans) :: this
  handle = atlas__Trans__handle (this%c_ptr())
end function

function truncation( this )
  use atlas_trans_c_binding
  integer :: truncation
  class(atlas_Trans) :: this
  truncation = atlas__Trans__truncation (this%c_ptr())
end function

function nb_spectral_coefficients( this )
  use atlas_trans_c_binding
  integer :: nb_spectral_coefficients
  class(atlas_Trans) :: this
  nb_spectral_coefficients = atlas__Trans__nspec2 (this%c_ptr())
end function

function nb_spectral_coefficients_global( this )
  use atlas_trans_c_binding
  integer :: nb_spectral_coefficients_global
  class(atlas_Trans) :: this
  nb_spectral_coefficients_global = atlas__Trans__nspec2g (this%c_ptr())
end function

function nb_gridpoints( this )
  use atlas_trans_c_binding
  integer :: nb_gridpoints
  class(atlas_Trans) :: this
  nb_gridpoints = atlas__Trans__ngptot (this%c_ptr())
end function

function nb_gridpoints_global( this )
  use atlas_trans_c_binding
  integer :: nb_gridpoints_global
  class(atlas_Trans) :: this
  nb_gridpoints_global = atlas__Trans__ngptotg (this%c_ptr())
end function

function grid( this )
  use atlas_trans_c_binding
  use atlas_grid_module
  class(atlas_Trans) :: this
  type(atlas_StructuredGrid) :: grid
  grid = atlas_StructuredGrid( atlas__Trans__grid(this%c_ptr()) )
  call grid%return()
end function

subroutine dirtrans_fieldset_nodes(this, gp, gpfields, sp, spfields, parameters)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_fieldset_module, only: atlas_FieldSet
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: gp
  class(atlas_FieldSet), intent(in)  :: gpfields
  class(atlas_FunctionSpace), intent(in)  :: sp
  class(atlas_FieldSet), intent(inout) :: spfields
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_fieldset_nodes( this%c_ptr(),     &
    &                          gp%c_ptr(), &
    &                          gpfields%c_ptr(), &
    &                          sp%c_ptr(), &
    &                          spfields%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine dirtrans_fieldset_nodes

subroutine dirtrans_fieldset(this, gpfields, spfields, parameters)
  use atlas_trans_c_binding
  use atlas_fieldset_module, only: atlas_FieldSet
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: gpfields
  class(atlas_FieldSet), intent(inout) :: spfields
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_fieldset( this%c_ptr(),     &
    &                          gpfields%c_ptr(), &
    &                          spfields%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine dirtrans_fieldset

subroutine invtrans_fieldset_nodes(this, sp, spfields, gp, gpfields, parameters)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_fieldset_module, only: atlas_FieldSet
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: sp
  class(atlas_FieldSet), intent(in)  :: spfields
  class(atlas_FunctionSpace), intent(in) :: gp
  class(atlas_FieldSet), intent(inout) :: gpfields
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_fieldset_nodes( this%c_ptr(),     &
    &                          sp%c_ptr(), &
    &                          spfields%c_ptr(), &
    &                          gp%c_ptr(), &
    &                          gpfields%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine invtrans_fieldset_nodes

subroutine invtrans_fieldset(this, spfields, gpfields, parameters)
  use atlas_trans_c_binding
  use atlas_fieldset_module, only: atlas_FieldSet
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: spfields
  class(atlas_FieldSet), intent(inout) :: gpfields
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_fieldset( this%c_ptr(),     &
    &                          spfields%c_ptr(), &
    &                          gpfields%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine invtrans_fieldset


subroutine dirtrans_field_nodes(this, gp, gpfield, sp, spfield, parameters)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: gp
  class(atlas_Field), intent(in)  :: gpfield
  class(atlas_FunctionSpace), intent(in) :: sp
  class(atlas_Field), intent(inout) :: spfield
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_field_nodes( this%c_ptr(), &
    &                          gp%c_ptr(), &
    &                          gpfield%c_ptr(), &
    &                          sp%c_ptr(), &
    &                          spfield%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine dirtrans_field_nodes

subroutine dirtrans_field(this, gpfield, spfield, parameters)
  use atlas_trans_c_binding
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: gpfield
  class(atlas_Field), intent(inout) :: spfield
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_field( this%c_ptr(), &
    &                          gpfield%c_ptr(), &
    &                          spfield%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine dirtrans_field

subroutine dirtrans_wind2vordiv_field(this, gp, gpwind, sp, spvor, spdiv, parameters)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: gp
  type(atlas_Field), intent(in)  :: gpwind
  class(atlas_FunctionSpace), intent(in) :: sp
  type(atlas_Field), intent(inout) :: spvor
  type(atlas_Field), intent(inout) :: spdiv
  type(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_wind2vordiv_field_nodes( this%c_ptr(), &
    &                          gp%c_ptr(), &
    &                          gpwind%c_ptr(), &
    &                          sp%c_ptr(), &
    &                          spvor%c_ptr(), &
    &                          spdiv%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine dirtrans_wind2vordiv_field

subroutine invtrans_field_nodes(this, sp, spfield, gp, gpfield, parameters)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: sp
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_FunctionSpace), intent(in)  :: gp
  class(atlas_Field), intent(inout) :: gpfield
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_field_nodes( this%c_ptr(), &
    &                          sp%c_ptr(), &
    &                          spfield%c_ptr(), &
    &                          gp%c_ptr(), &
    &                          gpfield%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine invtrans_field_nodes

subroutine invtrans_field(this, spfield, gpfield, parameters)
  use atlas_trans_c_binding
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_Field), intent(inout) :: gpfield
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_field( this%c_ptr(), &
    &                          spfield%c_ptr(), &
    &                          gpfield%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine invtrans_field


subroutine invtrans_vordiv2wind_field(this, sp, spvor, spdiv, gp, gpwind, parameters)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: sp
  class(atlas_Field), intent(in)  :: spvor
  class(atlas_Field), intent(in)  :: spdiv
  class(atlas_FunctionSpace), intent(in)  :: gp
  class(atlas_Field), intent(inout) :: gpwind
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    call p%reset_c_ptr( parameters%c_ptr() )
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_vordiv2wind_field_nodes( this%c_ptr(), &
    &                          sp%c_ptr(), &
    &                          spvor%c_ptr(), &
    &                          spdiv%c_ptr(), &
    &                          gp%c_ptr(), &
    &                          gpwind%c_ptr(), &
    &                          p%c_ptr() )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine invtrans_vordiv2wind_field


subroutine invtrans_grad_field_nodes(this, sp, spfield, gp, gpfield)
  use atlas_trans_c_binding
  use atlas_functionspace_module, only: atlas_Functionspace
  use atlas_field_module, only: atlas_Field
  class(atlas_Trans), intent(in) :: this
  class(atlas_FunctionSpace), intent(in)  :: sp
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_FunctionSpace), intent(in)  :: gp
  class(atlas_Field), intent(inout) :: gpfield
#ifdef ATLAS_HAVE_TRANS

  call atlas__Trans__invtrans_grad_field_nodes( this%c_ptr(), &
    &                          sp%c_ptr(), &
    &                          spfield%c_ptr(), &
    &                          gp%c_ptr(), &
    &                          gpfield%c_ptr() )
#else
  THROW_ERROR
#endif
end subroutine invtrans_grad_field_nodes



subroutine gathspec_r1(this, local, global)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: local(:)
  real(c_double), intent(inout) :: global(:)
#ifdef ATLAS_HAVE_TRANS
  call atlas__Trans__gathspec(this%c_ptr(), 1, (/1/), local, global )
#else
  THROW_ERROR
#endif
end subroutine gathspec_r1

subroutine gathspec_r2(this, local, global)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_array_module, only: array_view1d
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: local(:,:)
  real(c_double), intent(inout) :: global(:,:)
#ifdef ATLAS_HAVE_TRANS
  real(c_double), pointer :: local_view(:), global_view(:)
  integer :: destination(size(local,1))
  destination(:) = 1
  local_view => array_view1d(local)
  global_view => array_view1d(global)
  call atlas__Trans__gathspec(this%c_ptr(), size(local,1), destination, local_view, global_view )
#else
  THROW_ERROR
#endif
end subroutine gathspec_r2


subroutine specnorm_r1_scalar(this, spectra, norm, rank)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: spectra(:)
  real(c_double), intent(out) :: norm
  integer, optional :: rank
#ifdef ATLAS_HAVE_TRANS
  integer :: rank_opt
  real(c_double) :: norms(1)
  rank_opt = 0
  if( present(rank) ) rank_opt = rank
  call atlas__Trans__specnorm(this%c_ptr(), 1, spectra, norms, rank_opt )
  norm = norms(1)
#else
  norm=0
  THROW_ERROR
#endif
end subroutine

subroutine specnorm_r2(this, spectra, norm, rank)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_array_module, only: array_view1d
  use atlas_trans_c_binding
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: spectra(:,:)
  real(c_double), intent(inout) :: norm(:)
  integer, optional :: rank
#ifdef ATLAS_HAVE_TRANS
  integer :: rank_opt
  real(c_double), pointer :: spectra_view(:)
  spectra_view => array_view1d(spectra)
  call atlas__Trans__specnorm(this%c_ptr(), size(spectra,1), spectra_view, norm, rank_opt )
#else
  THROW_ERROR
#endif
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_Trans_module
