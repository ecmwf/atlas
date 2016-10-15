
#include "atlas/atlas_f.h"

module atlas_Trans_module


use fckit_c_interop_module, only: c_str
use fckit_array_module, only: array_view1d
use fckit_object_module, only: fckit_object
use fckit_refcounted_module, only: fckit_refcounted
use atlas_Grid_module, only: atlas_Grid
use atlas_functionspace_module, only: atlas_Functionspace
use atlas_field_module, only: atlas_Field
use atlas_fieldset_module, only: atlas_FieldSet
use atlas_Error_module, only: atlas_code_location, atlas_throw_usererror

implicit none

private :: c_str, array_view1d
private :: fckit_refcounted
private :: fckit_object
private :: atlas_Grid
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_FunctionSpace
private :: atlas_code_location

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
  procedure :: handle       => atlas_Trans__handle
  procedure :: nproc        => atlas_Trans__nproc
  procedure :: myproc       => atlas_Trans__myproc
  procedure :: ndgl         => atlas_Trans__ndgl
  procedure :: nsmax        => atlas_Trans__nsmax
  procedure :: ngptot       => atlas_Trans__ngptot
  procedure :: ngptotg      => atlas_Trans__ngptotg
  procedure :: ngptotmx     => atlas_Trans__ngptotmx
  procedure :: nspec        => atlas_Trans__nspec
  procedure :: nspec2       => atlas_Trans__nspec2
  procedure :: nspec2g      => atlas_Trans__nspec2g
  procedure :: nspec2mx     => atlas_Trans__nspec2mx
  procedure :: n_regions_NS => atlas_Trans__n_regions_NS
  procedure :: n_regions_EW => atlas_Trans__n_regions_EW
  procedure :: nloen        => atlas_Trans__nloen
  procedure :: n_regions    => atlas_Trans__n_regions
  procedure :: nfrstlat     => atlas_Trans__nfrstlat
  procedure :: nlstlat      => atlas_Trans__nlstlat
  procedure :: nptrfrstlat  => atlas_Trans__nptrfrstlat
  procedure :: nsta         => atlas_Trans__nsta
  procedure :: nonl         => atlas_Trans__nonl
  procedure :: nmyms        => atlas_Trans__nmyms
  procedure :: nasm0        => atlas_Trans__nasm0
  procedure :: nump         => atlas_Trans__nump
  procedure :: nvalue       => atlas_Trans__nvalue

  procedure, private :: dirtrans_field_nodes => atlas_Trans__dirtrans_field_nodes
  procedure, private :: dirtrans_field => atlas_Trans__dirtrans_field
  procedure, private :: dirtrans_fieldset_nodes => atlas_Trans__dirtrans_fieldset_nodes
  procedure, private :: dirtrans_fieldset => atlas_Trans__dirtrans_fieldset
  procedure, public :: dirtrans_wind2vordiv => atlas_Trans__dirtrans_wind2vordiv_field
  generic, public :: dirtrans => &
    & dirtrans_field, &
    & dirtrans_fieldset, &
    & dirtrans_fieldset_nodes, &
    & dirtrans_field_nodes

  procedure, private :: invtrans_field_nodes => atlas_Trans__invtrans_field_nodes
  procedure, private :: invtrans_field => atlas_Trans__invtrans_field
  procedure, private :: invtrans_fieldset_nodes => atlas_Trans__invtrans_fieldset_nodes
  procedure, private :: invtrans_fieldset => atlas_Trans__invtrans_fieldset
  procedure, public :: invtrans_vordiv2wind => atlas_Trans__invtrans_vordiv2wind_field
  generic, public :: invtrans => &
    & invtrans_field, &
    & invtrans_fieldset, &
    & invtrans_field_nodes, &
    & invtrans_fieldset_nodes

  procedure, private :: gathspec_r1 => atlas_Trans__gathspec_r1
  procedure, private :: gathspec_r2 => atlas_Trans__gathspec_r2
  generic, public :: gathspec => gathspec_r1, gathspec_r2

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
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  call atlas_throw_usererror("Cannot use atlas_Trans since atlas is compiled without" // &
    & "ENABLE_TRANS=ON",atlas_code_location(file,line))
end subroutine

function atlas_Trans__ctor( grid, nsmax ) result(trans)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only: c_null_ptr
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

function atlas_Trans__handle( this ) result(handle)
  use atlas_trans_c_binding
  integer :: handle
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  handle = atlas__Trans__handle (this%c_ptr())
#else
  THROW_ERROR
  handle = 0
#endif
end function


function atlas_Trans__nproc( this ) result(nproc)
  use atlas_trans_c_binding
  integer :: nproc
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nproc = atlas__Trans__nproc (this%c_ptr())
#else
  THROW_ERROR
  nproc = 0
#endif
end function


function atlas_Trans__myproc( this, proc0 ) result(myproc)
  use atlas_trans_c_binding
  integer :: myproc
  class(atlas_Trans) :: this
  integer, intent(in) :: proc0
#ifdef ATLAS_HAVE_TRANS
  myproc = atlas__Trans__myproc (this%c_ptr(),proc0)
#else
  myproc = 0
#endif
end function

function atlas_Trans__ndgl( this ) result(ndgl)
  use atlas_trans_c_binding
  integer :: ndgl
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ndgl = atlas__Trans__ndgl (this%c_ptr())
#else
  THROW_ERROR
  ndgl = 0
#endif
end function

function atlas_Trans__nsmax( this ) result(nsmax)
  use atlas_trans_c_binding
  integer :: nsmax
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nsmax = atlas__Trans__nsmax (this%c_ptr())
#else
  nsmax = 0
#endif
end function

function atlas_Trans__ngptot( this ) result(ngptot)
  use atlas_trans_c_binding
  integer :: ngptot
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ngptot = atlas__Trans__ngptot (this%c_ptr())
#else
  THROW_ERROR
  ngptot = 0
#endif
end function

function atlas_Trans__ngptotg( this ) result(ngptotg)
  use atlas_trans_c_binding
  integer :: ngptotg
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ngptotg = atlas__Trans__ngptotg (this%c_ptr())
#else
  THROW_ERROR
  ngptotg = 0
#endif
end function

function atlas_Trans__ngptotmx( this ) result(ngptotmx)
  use atlas_trans_c_binding
  integer :: ngptotmx
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ngptotmx = atlas__Trans__ngptotmx (this%c_ptr())
#else
  THROW_ERROR
  ngptotmx = 0
#endif
end function


function atlas_Trans__nspec( this ) result(nspec)
  use atlas_trans_c_binding
  integer :: nspec
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec = atlas__Trans__nspec (this%c_ptr())
#else
  nspec = 0
#endif
end function

function atlas_Trans__nspec2( this ) result(nspec2)
  use atlas_trans_c_binding
  integer :: nspec2
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec2 = atlas__Trans__nspec2 (this%c_ptr())
#else
  nspec2 = 0
#endif
end function

function atlas_Trans__nspec2g( this ) result(nspec2g)
  use atlas_trans_c_binding
  integer :: nspec2g
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec2g = atlas__Trans__nspec2g (this%c_ptr())
#else
  nspec2g = 0
#endif
end function

function atlas_Trans__nspec2mx( this ) result(nspec2mx)
  use atlas_trans_c_binding
  integer :: nspec2mx
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec2mx = atlas__Trans__nspec2mx (this%c_ptr())
#else
  nspec2mx = 0
#endif
end function


function atlas_Trans__n_regions_NS( this ) result(n_regions_NS)
  use atlas_trans_c_binding
  integer :: n_regions_NS
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  n_regions_NS = atlas__Trans__n_regions_NS (this%c_ptr())
#else
  THROW_ERROR
  n_regions_NS = 0
#endif
end function


function atlas_Trans__n_regions_EW( this ) result(n_regions_EW)
  use atlas_trans_c_binding
  integer :: n_regions_EW
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  n_regions_EW = atlas__Trans__n_regions_EW (this%c_ptr())
#else
  THROW_ERROR
  n_regions_EW = 0
#endif
end function

function atlas_Trans__nump( this ) result(nump)
  use atlas_trans_c_binding
  integer :: nump
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nump = atlas__Trans__nump (this%c_ptr())
#else
  THROW_ERROR
  nump = 0
#endif
end function

function atlas_Trans__nloen(this) result(nloen)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nloen(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nloen_c_ptr
  integer(c_int) :: size
  nloen_c_ptr =  atlas__Trans__nloen(this%c_ptr(), size)
  call c_f_pointer ( nloen_c_ptr , nloen , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nloen) ) then
    allocate(nloen(0) )
  endif
#endif
end function atlas_Trans__nloen

function atlas_Trans__n_regions(this) result(n_regions)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: n_regions(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: n_regions_c_ptr
  integer(c_int) :: size
  n_regions_c_ptr =  atlas__Trans__n_regions(this%c_ptr(), size)
  call c_f_pointer ( n_regions_c_ptr , n_regions , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(n_regions) ) then
    allocate(n_regions(0) )
  endif
#endif
end function atlas_Trans__n_regions


function atlas_Trans__nfrstlat(this) result(nfrstlat)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nfrstlat(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nfrstlat_c_ptr
  integer(c_int) :: size
  nfrstlat_c_ptr =  atlas__Trans__nfrstlat(this%c_ptr(), size)
  call c_f_pointer ( nfrstlat_c_ptr , nfrstlat , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nfrstlat) ) then
    allocate(nfrstlat(0) )
  endif
#endif
end function atlas_Trans__nfrstlat

function atlas_Trans__nlstlat(this) result(nlstlat)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_f_pointer
  integer(c_int), pointer :: nlstlat(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nlstlat_c_ptr
  integer(c_int) :: size
  nlstlat_c_ptr =  atlas__Trans__nlstlat(this%c_ptr(), size)
  call c_f_pointer ( nlstlat_c_ptr , nlstlat , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nlstlat) ) then
    allocate(nlstlat(0) )
  endif
#endif
end function atlas_Trans__nlstlat


function atlas_Trans__nptrfrstlat(this) result(nptrfrstlat)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nptrfrstlat(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nptrfrstlat_c_ptr
  integer(c_int) :: size
  nptrfrstlat_c_ptr =  atlas__Trans__nptrfrstlat(this%c_ptr(), size)
  call c_f_pointer ( nptrfrstlat_c_ptr , nptrfrstlat , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nptrfrstlat) ) then
    allocate(nptrfrstlat(0) )
  endif
#endif
end function atlas_Trans__nptrfrstlat


function atlas_Trans__nsta(this) result(nsta)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nsta(:,:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nsta_c_ptr
  integer(c_int) :: sizef1
  integer(c_int) :: sizef2
  nsta_c_ptr =  atlas__Trans__nsta(this%c_ptr(), sizef2, sizef1)
  call c_f_pointer ( nsta_c_ptr , nsta , (/sizef1,sizef2/) )
#else
  THROW_ERROR
  if ( .not. associated(nsta) ) then
    allocate(nsta(0,0) )
  endif
#endif
end function atlas_Trans__nsta

function atlas_Trans__nonl(this) result(nonl)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nonl(:,:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nonl_c_ptr
  integer(c_int) :: sizef1
  integer(c_int) :: sizef2
  nonl_c_ptr =  atlas__Trans__nonl(this%c_ptr(), sizef2, sizef1)
  call c_f_pointer ( nonl_c_ptr , nonl , (/sizef1,sizef2/) )
#else
  THROW_ERROR
  if ( .not. associated(nonl) ) then
    allocate(nonl(0,0) )
  endif
#endif
end function atlas_Trans__nonl


function atlas_Trans__nmyms(this) result(nmyms)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nmyms(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nmyms_c_ptr
  integer(c_int) :: size
  nmyms_c_ptr =  atlas__Trans__nptrfrstlat(this%c_ptr(), size)
  call c_f_pointer ( nmyms_c_ptr , nmyms , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nmyms) ) then
    allocate(nmyms(0) )
  endif
#endif
end function atlas_Trans__nmyms

function atlas_Trans__nasm0(this) result(nasm0)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nasm0(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nasm0_c_ptr
  integer(c_int) :: size
  nasm0_c_ptr =  atlas__Trans__nasm0(this%c_ptr(), size)
  call c_f_pointer ( nasm0_c_ptr , nasm0 , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nasm0) ) then
    allocate(nasm0(0) )
  endif
#endif
end function atlas_Trans__nasm0

function atlas_Trans__nvalue(this) result(nvalue)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_int, c_f_pointer, c_ptr
  integer(c_int), pointer :: nvalue(:)
  class(atlas_Trans), intent(in) :: this
#ifdef ATLAS_HAVE_TRANS
  type(c_ptr) :: nvalue_c_ptr
  integer(c_int) :: size
  nvalue_c_ptr =  atlas__Trans__nvalue(this%c_ptr(), size)
  call c_f_pointer ( nvalue_c_ptr , nvalue , (/size/) )
#else
  THROW_ERROR
  if ( .not. associated(nvalue) ) then
    allocate(nvalue(0) )
  endif
#endif
end function atlas_Trans__nvalue

subroutine atlas_Trans__dirtrans_fieldset_nodes(this, gp, gpfields, sp, spfields, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__dirtrans_fieldset_nodes

subroutine atlas_Trans__dirtrans_fieldset(this, gpfields, spfields, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__dirtrans_fieldset

subroutine atlas_Trans__invtrans_fieldset_nodes(this, sp, spfields, gp, gpfields, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__invtrans_fieldset_nodes

subroutine atlas_Trans__invtrans_fieldset(this, spfields, gpfields, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__invtrans_fieldset


subroutine atlas_Trans__dirtrans_field_nodes(this, gp, gpfield, sp, spfield, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__dirtrans_field_nodes

subroutine atlas_Trans__dirtrans_field(this, gpfield, spfield, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__dirtrans_field

subroutine atlas_Trans__dirtrans_wind2vordiv_field(this, gp, gpwind, sp, spvor, spdiv, parameters)
  use atlas_trans_c_binding
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

end subroutine atlas_Trans__dirtrans_wind2vordiv_field

subroutine atlas_Trans__invtrans_field_nodes(this, sp, spfield, gp, gpfield, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__invtrans_field_nodes

subroutine atlas_Trans__invtrans_field(this, spfield, gpfield, parameters)
  use atlas_trans_c_binding
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
end subroutine atlas_Trans__invtrans_field


subroutine atlas_Trans__invtrans_vordiv2wind_field(this, sp, spvor, spdiv, gp, gpwind, parameters)
  use atlas_trans_c_binding
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

end subroutine atlas_Trans__invtrans_vordiv2wind_field


subroutine atlas_Trans__gathspec_r1(this, local, global)
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
end subroutine atlas_Trans__gathspec_r1

subroutine atlas_Trans__gathspec_r2(this, local, global)
  use atlas_trans_c_binding
  use, intrinsic :: iso_c_binding, only : c_double
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
end subroutine atlas_Trans__gathspec_r2

! ----------------------------------------------------------------------------------------

end module atlas_Trans_module
