! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! Trans routines

#ifdef ATLAS_HAVE_TRANS
#define USE_ATLAS_TRANS_C_BINDING   use atlas_trans_c_binding
#else
#define USE_ATLAS_TRANS_C_BINDING
#define THROW_ERROR call atlas_throw_usererror("Cannot use atlas_Trans since atlas is compiled without ENABLE_TRANS=ON",atlas_code_location(__FILE__,__LINE__))
#endif

function atlas_Trans__ctor( grid, nsmax ) result(trans)
  USE_ATLAS_TRANS_C_BINDING
  type(atlas_Trans) :: trans
  type(atlas_ReducedGrid), intent(in) :: grid
  integer, intent(in), optional :: nsmax
#ifdef ATLAS_HAVE_TRANS
  if( present(nsmax) ) then
    trans%cpp_object_ptr = atlas__Trans__new( grid%cpp_object_ptr, nsmax )
  else
    trans%cpp_object_ptr = atlas__Trans__new( grid%cpp_object_ptr, 0 )
  endif
#else
  ! IGNORE
#endif
end function atlas_Trans__ctor

function atlas_TransParameters__ctor() result(params)
  USE_ATLAS_TRANS_C_BINDING
  type(atlas_TransParameters) :: params
#ifdef ATLAS_HAVE_TRANS
  params%cpp_object_ptr = atlas__TransParameters__new()
#else
  ! IGNORE
#endif
end function atlas_TransParameters__ctor

subroutine atlas_Trans__delete( trans )
  USE_ATLAS_TRANS_C_BINDING
  type(atlas_Trans) :: trans
#ifdef ATLAS_HAVE_TRANS
  call atlas__Trans__delete(trans%cpp_object_ptr);
#else
  ! IGNORE
#endif
end subroutine

subroutine atlas_TransParameters__delete( parameters )
  USE_ATLAS_TRANS_C_BINDING
  type(atlas_TransParameters) :: parameters
#ifdef ATLAS_HAVE_TRANS
  call atlas__TransParameters__delete(parameters%cpp_object_ptr);
#else
  ! IGNORE
#endif
end subroutine

function atlas_Trans__handle( this ) result(handle)
  USE_ATLAS_TRANS_C_BINDING
  integer :: handle
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  handle = atlas__Trans__handle (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function


function atlas_Trans__nproc( this ) result(nproc)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nproc
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nproc = atlas__Trans__nproc (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function


function atlas_Trans__myproc( this, proc0 ) result(myproc)
  USE_ATLAS_TRANS_C_BINDING
  integer :: myproc
  class(atlas_Trans) :: this
  integer, intent(in) :: proc0
#ifdef ATLAS_HAVE_TRANS
  myproc = atlas__Trans__myproc (this%cpp_object_ptr,proc0)
#else
  THROW_ERROR
#endif
end function

function atlas_Trans__ndgl( this ) result(ndgl)
  USE_ATLAS_TRANS_C_BINDING
  integer :: ndgl
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ndgl = atlas__Trans__ndgl (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function

function atlas_Trans__nsmax( this ) result(nsmax)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nsmax
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nsmax = atlas__Trans__nsmax (this%cpp_object_ptr)
#else
  nsmax = 0
#endif
end function

function atlas_Trans__ngptot( this ) result(ngptot)
  USE_ATLAS_TRANS_C_BINDING
  integer :: ngptot
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ngptot = atlas__Trans__ngptot (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function

function atlas_Trans__ngptotg( this ) result(ngptotg)
  USE_ATLAS_TRANS_C_BINDING
  integer :: ngptotg
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ngptotg = atlas__Trans__ngptotg (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function

function atlas_Trans__ngptotmx( this ) result(ngptotmx)
  USE_ATLAS_TRANS_C_BINDING
  integer :: ngptotmx
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  ngptotmx = atlas__Trans__ngptotmx (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function


function atlas_Trans__nspec( this ) result(nspec)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nspec
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec = atlas__Trans__nspec (this%cpp_object_ptr)
#else
  nspec = 0
#endif
end function

function atlas_Trans__nspec2( this ) result(nspec2)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nspec2
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec2 = atlas__Trans__nspec2 (this%cpp_object_ptr)
#else
  nspec2 = 0
#endif
end function

function atlas_Trans__nspec2g( this ) result(nspec2g)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nspec2g
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec2g = atlas__Trans__nspec2g (this%cpp_object_ptr)
#else
  nspec2g = 0
#endif
end function

function atlas_Trans__nspec2mx( this ) result(nspec2mx)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nspec2mx
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nspec2mx = atlas__Trans__nspec2mx (this%cpp_object_ptr)
#else
  nspec2mx = 0
#endif
end function


function atlas_Trans__n_regions_NS( this ) result(n_regions_NS)
  USE_ATLAS_TRANS_C_BINDING
  integer :: n_regions_NS
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  n_regions_NS = atlas__Trans__n_regions_NS (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function


function atlas_Trans__n_regions_EW( this ) result(n_regions_EW)
  USE_ATLAS_TRANS_C_BINDING
  integer :: n_regions_EW
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  n_regions_EW = atlas__Trans__n_regions_EW (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function

function atlas_Trans__nump( this ) result(nump)
  USE_ATLAS_TRANS_C_BINDING
  integer :: nump
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  nump = atlas__Trans__nump (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function

function atlas_Trans__nloen(this) result(nloen)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nloen(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nloen_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nloen_c_ptr =  atlas__Trans__nloen(this%cpp_object_ptr, size)
  call C_F_POINTER ( nloen_c_ptr , nloen , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nloen

function atlas_Trans__n_regions(this) result(n_regions)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: n_regions(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: n_regions_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  n_regions_c_ptr =  atlas__Trans__n_regions(this%cpp_object_ptr, size)
  call C_F_POINTER ( n_regions_c_ptr , n_regions , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__n_regions


function atlas_Trans__nfrstlat(this) result(nfrstlat)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nfrstlat(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nfrstlat_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nfrstlat_c_ptr =  atlas__Trans__nfrstlat(this%cpp_object_ptr, size)
  call C_F_POINTER ( nfrstlat_c_ptr , nfrstlat , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nfrstlat

function atlas_Trans__nlstlat(this) result(nlstlat)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nlstlat(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nlstlat_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nlstlat_c_ptr =  atlas__Trans__nlstlat(this%cpp_object_ptr, size)
  call C_F_POINTER ( nlstlat_c_ptr , nlstlat , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nlstlat


function atlas_Trans__nptrfrstlat(this) result(nptrfrstlat)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nptrfrstlat(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nptrfrstlat_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nptrfrstlat_c_ptr =  atlas__Trans__nptrfrstlat(this%cpp_object_ptr, size)
  call C_F_POINTER ( nptrfrstlat_c_ptr , nptrfrstlat , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nptrfrstlat


function atlas_Trans__nsta(this) result(nsta)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nsta(:,:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nsta_c_ptr
  integer(c_int) :: sizef1
  integer(c_int) :: sizef2
#ifdef ATLAS_HAVE_TRANS
  nsta_c_ptr =  atlas__Trans__nsta(this%cpp_object_ptr, sizef2, sizef1)
  call C_F_POINTER ( nsta_c_ptr , nsta , (/sizef1,sizef2/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nsta

function atlas_Trans__nonl(this) result(nonl)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nonl(:,:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nonl_c_ptr
  integer(c_int) :: sizef1
  integer(c_int) :: sizef2
#ifdef ATLAS_HAVE_TRANS
  nonl_c_ptr =  atlas__Trans__nonl(this%cpp_object_ptr, sizef2, sizef1)
  call C_F_POINTER ( nonl_c_ptr , nonl , (/sizef1,sizef2/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nonl


function atlas_Trans__nmyms(this) result(nmyms)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nmyms(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nmyms_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nmyms_c_ptr =  atlas__Trans__nptrfrstlat(this%cpp_object_ptr, size)
  call C_F_POINTER ( nmyms_c_ptr , nmyms , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nmyms

function atlas_Trans__nasm0(this) result(nasm0)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nasm0(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nasm0_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nasm0_c_ptr =  atlas__Trans__nasm0(this%cpp_object_ptr, size)
  call C_F_POINTER ( nasm0_c_ptr , nasm0 , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nasm0

function atlas_Trans__nvalue(this) result(nvalue)
  USE_ATLAS_TRANS_C_BINDING
  integer(c_int), pointer :: nvalue(:)
  class(atlas_Trans), intent(in) :: this
  type(c_ptr) :: nvalue_c_ptr
  integer(c_int) :: size
#ifdef ATLAS_HAVE_TRANS
  nvalue_c_ptr =  atlas__Trans__nvalue(this%cpp_object_ptr, size)
  call C_F_POINTER ( nvalue_c_ptr , nvalue , (/size/) )
#else
  THROW_ERROR
#endif
end function atlas_Trans__nvalue

subroutine atlas_Trans__dirtrans_fieldset(this, gpfields, spfields, parameters)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: gpfields
  class(atlas_FieldSet), intent(inout) :: spfields
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    p%cpp_object_ptr = parameters%cpp_object_ptr
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_fieldset( this%cpp_object_ptr,     &
    &                          gpfields%cpp_object_ptr, &
    &                          spfields%cpp_object_ptr, &
    &                          p%cpp_object_ptr )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif
end subroutine atlas_Trans__dirtrans_fieldset

subroutine atlas_Trans__invtrans_fieldset(this, spfields, gpfields, parameters)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  class(atlas_FieldSet), intent(in)  :: spfields
  class(atlas_FieldSet), intent(inout) :: gpfields
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    p%cpp_object_ptr = parameters%cpp_object_ptr
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_fieldset( this%cpp_object_ptr,     &
    &                          spfields%cpp_object_ptr, &
    &                          gpfields%cpp_object_ptr, &
    &                          p%cpp_object_ptr )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine atlas_Trans__invtrans_fieldset

subroutine atlas_Trans__dirtrans_field(this, gpfield, spfield, parameters)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: gpfield
  class(atlas_Field), intent(inout) :: spfield
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    p%cpp_object_ptr = parameters%cpp_object_ptr
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_field( this%cpp_object_ptr, &
    &                          gpfield%cpp_object_ptr, &
    &                          spfield%cpp_object_ptr, &
    &                          p%cpp_object_ptr )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine atlas_Trans__dirtrans_field

subroutine atlas_Trans__dirtrans_wind2vordiv_field(this, gpwind, spvor, spdiv, parameters)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: gpwind
  class(atlas_Field), intent(inout) :: spvor
  class(atlas_Field), intent(inout) :: spdiv
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    p%cpp_object_ptr = parameters%cpp_object_ptr
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__dirtrans_wind2vordiv_field( this%cpp_object_ptr, &
    &                          gpwind%cpp_object_ptr, &
    &                          spvor%cpp_object_ptr, &
    &                          spdiv%cpp_object_ptr, &
    &                          p%cpp_object_ptr )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine atlas_Trans__dirtrans_wind2vordiv_field

subroutine atlas_Trans__invtrans_field(this, spfield, gpfield, parameters)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spfield
  class(atlas_Field), intent(inout) :: gpfield
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    p%cpp_object_ptr = parameters%cpp_object_ptr
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_field( this%cpp_object_ptr, &
    &                          spfield%cpp_object_ptr, &
    &                          gpfield%cpp_object_ptr, &
    &                          p%cpp_object_ptr )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine atlas_Trans__invtrans_field

subroutine atlas_Trans__invtrans_vordiv2wind_field(this, spvor, spdiv, gpwind, parameters)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  class(atlas_Field), intent(in)  :: spvor
  class(atlas_Field), intent(in)  :: spdiv
  class(atlas_Field), intent(inout) :: gpwind
  class(atlas_TransParameters), intent(in), optional  :: parameters
#ifdef ATLAS_HAVE_TRANS
  type(atlas_TransParameters) :: p

  if( present(parameters) ) then
    p%cpp_object_ptr = parameters%cpp_object_ptr
  else
    p = atlas_TransParameters()
  endif

  call atlas__Trans__invtrans_vordiv2wind_field( this%cpp_object_ptr, &
    &                          spvor%cpp_object_ptr, &
    &                          spdiv%cpp_object_ptr, &
    &                          gpwind%cpp_object_ptr, &
    &                          p%cpp_object_ptr )

  if( .not. present(parameters) ) then
    call atlas_TransParameters__delete(p)
  endif
#else
  THROW_ERROR
#endif

end subroutine atlas_Trans__invtrans_vordiv2wind_field


subroutine atlas_Trans__gathspec_r1(this, local, global)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: local(:)
  real(c_double), intent(out) :: global(:)
#ifdef ATLAS_HAVE_TRANS
  call atlas__Trans__gathspec(this%cpp_object_ptr, 1, (/1/), local, global )
#endif
end subroutine atlas_Trans__gathspec_r1

subroutine atlas_Trans__gathspec_r2(this, local, global)
  USE_ATLAS_TRANS_C_BINDING
  class(atlas_Trans), intent(in) :: this
  real(c_double), intent(in) :: local(:,:)
  real(c_double), intent(out) :: global(:,:)
#ifdef ATLAS_HAVE_TRANS
  integer :: destination(size(local,2))
  destination(:) = 1
  call atlas__Trans__gathspec(this%cpp_object_ptr, size(local,2), destination, view1d(local), view1d(global) )
#endif
end subroutine atlas_Trans__gathspec_r2
