
module atlas_Grid_module


use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_Grid
public :: atlas_StructuredGrid
public :: atlas_GaussianGrid
public :: atlas_ReducedGaussianGrid
public :: atlas_RegularGaussianGrid
public :: atlas_RegularLonLatGrid

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_Grid

! Purpose :
! -------
!   *atlas_Grid* : Object Grid specifications for Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: npts => atlas_Grid__size
  procedure, public :: delete => atlas_Grid__delete
  procedure, public :: copy => atlas_Grid__copy
END TYPE atlas_Grid

!------------------------------------------------------------------------------

TYPE, extends(atlas_Grid) :: atlas_StructuredGrid

! Purpose :
! -------
!   *atlas_StructuredGrid* : Object Grid specifications for Reduced Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: ny        => Structured__ny
  procedure, private   :: nx_int32 => Structured__nx_int32
  procedure, private   :: nx_int64 => Structured__nx_int64
  generic   :: nx => nx_int32, nx_int64
  procedure :: nx_array  => Structured__nx_array
  procedure :: nxmin     => Structured__nxmin
  procedure :: nxmax     => Structured__nxmax
  procedure :: y_array   => Structured__y_array
  procedure, private :: x_32  => Structured__x_32
  procedure, private :: x_64  => Structured__x_64
  generic  :: x         => x_32, x_64
  procedure, private :: y_32         => Structured__y_32
  procedure, private :: y_64         => Structured__y_64
  generic :: y         => y_32, y_64
  procedure, private :: xy_32        => Structured__xy_32
  procedure, private :: xy_64        => Structured__xy_64
  generic :: xy        => xy_32, xy_64
  procedure, private :: lonlat_32    => Structured__lonlat_32
  procedure, private :: lonlat_64    => Structured__lonlat_64
  generic :: lonlat    => lonlat_32, lonlat_64
  procedure :: reduced   => Structured__reduced
END TYPE atlas_StructuredGrid

interface atlas_StructuredGrid
  module procedure atlas_StructuredGrid__ctor_id
  module procedure atlas_StructuredGrid__ctor_config
end interface

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_GaussianGrid

! Purpose :
! -------
!   *atlas_ReducedGaussianGrid* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
procedure :: N         => Gaussian__N
END TYPE atlas_GaussianGrid

interface atlas_GaussianGrid
  module procedure atlas_GaussianGrid__ctor_id
end interface

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_ReducedGaussianGrid

! Purpose :
! -------
!   *atlas_ReducedGaussianGrid* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
procedure :: N         => ReducedGaussian__N
END TYPE atlas_ReducedGaussianGrid

interface atlas_ReducedGaussianGrid
  module procedure atlas_ReducedGaussianGrid__ctor_N_int32_pl_int32
  module procedure atlas_ReducedGaussianGrid__ctor_N_int32_pl_int64
  module procedure atlas_ReducedGaussianGrid__ctor_N_int64_pl_int32
  module procedure atlas_ReducedGaussianGrid__ctor_N_int64_pl_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_RegularGaussianGrid

! Purpose :
! -------
!   *atlas_RegularGaussianGrid* : Object Grid specifications for Regular Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
procedure :: N         => RegularGaussian__N
END TYPE atlas_RegularGaussianGrid

interface atlas_RegularGaussianGrid
  module procedure atlas_RegularGaussianGrid__ctor_int32
  module procedure atlas_RegularGaussianGrid__ctor_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_RegularLonLatGrid

! Purpose :
! -------
!   *atlas_RegularLonLatGrid* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_RegularLonLatGrid

interface atlas_RegularLonLatGrid
  module procedure atlas_grid_RegularLonLat__ctor_int32
  module procedure atlas_grid_RegularLonLat__ctor_int64
end interface

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

interface c_idx
  module procedure :: c_idx_32
  module procedure :: c_idx_64
end interface

!------------------------------------------------------------------------------
!========================================================
contains
!========================================================

pure function c_idx_32(f_idx) result(c_idx)
    use, intrinsic :: iso_c_binding, only : c_long
    integer(c_long) :: c_idx
    integer(c_long), intent(in) :: f_idx
    c_idx = f_idx - 1_c_long
end function

pure function c_idx_64(f_idx) result(c_idx)
    use, intrinsic :: iso_c_binding, only : c_long, c_int
    integer(c_long) :: c_idx
    integer(c_int), intent(in) :: f_idx
    c_idx = f_idx - 1_c_long
end function


! -----------------------------------------------------------------------------
! Constructors

function atlas_StructuredGrid__ctor_id(identifier) result(grid)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_StructuredGrid) :: grid
  character(len=*), intent(in) :: identifier
  call grid%reset_c_ptr( atlas__grid__Structured(c_str(identifier)) )
end function

function atlas_StructuredGrid__ctor_config(config) result(grid)
  use atlas_grid_Structured_c_binding
  use atlas_Config_module, only: atlas_Config
  type(atlas_StructuredGrid) :: grid
  type(atlas_Config), intent(in) :: config
  call grid%reset_c_ptr( atlas__grid__Structured__config(config%c_ptr()) )
end function

!-----------------------------------------------------------------------------

function atlas_GaussianGrid__ctor_id(identifier) result(grid)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_GaussianGrid) :: grid
  character(len=*), intent(in) :: identifier
  call grid%reset_c_ptr( atlas__grid__Structured(c_str(identifier)) )
end function

! -----------------------------------------------------------------------------

function atlas_RegularGaussianGrid__ctor_int32(N) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularGaussianGrid) :: grid
  integer(c_int), intent(in) :: N
  call grid%reset_c_ptr( atlas__grid__regular__RegularGaussian(int(N,c_long)) )
end function

function atlas_RegularGaussianGrid__ctor_int64(N) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularGaussianGrid) :: grid
  integer(c_long), intent(in) :: N
  call grid%reset_c_ptr( atlas__grid__regular__RegularGaussian(int(N,c_long)) )
end function

!-----------------------------------------------------------------------------

function atlas_ReducedGaussianGrid__ctor_N_int32_pl_int32(N,pl) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: grid
  integer(c_int), intent(in) :: N
  integer(c_int), intent(in)  :: pl(:)
  call grid%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_int( int(N,c_long), pl) )
end function

function atlas_ReducedGaussianGrid__ctor_N_int32_pl_int64(N,pl) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: grid
  integer(c_int), intent(in) :: N
  integer(c_long), intent(in) :: pl(:)
  call grid%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_long(int(N,c_long),pl) )
end function

function atlas_ReducedGaussianGrid__ctor_N_int64_pl_int32(N,pl) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: grid
  integer(c_long), intent(in) :: N
  integer(c_int), intent(in)  :: pl(:)
  call grid%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_int(N, pl) )
end function

function atlas_ReducedGaussianGrid__ctor_N_int64_pl_int64(N,pl) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: grid
  integer(c_long), intent(in) :: N
  integer(c_long), intent(in) :: pl(:)
  call grid%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_long(N,pl) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_RegularLonLat__ctor_int32(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularLonLatGrid) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__RegularLonLat(int(nlon,c_long),int(nlat,c_long)) )
end function

function atlas_grid_RegularLonLat__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularLonLatGrid) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__RegularLonLat( nlon, nlat ) )
end function

! -----------------------------------------------------------------------------
! Structured members

function atlas_Grid__size(this) result(npts)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_Grid), intent(in) :: this
  integer(c_long) :: npts
  npts = atlas__grid__Structured__size(this%c_ptr())
end function

function Gaussian__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_GaussianGrid), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Gaussian__N(this%c_ptr())
end function

function ReducedGaussian__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_ReducedGaussianGrid), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Gaussian__N(this%c_ptr())
end function

function RegularGaussian__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_RegularGaussianGrid), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Gaussian__N(this%c_ptr())
end function

function Structured__ny(this) result(ny)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) :: ny
  ny = atlas__grid__Structured__ny(this%c_ptr())
end function


function Structured__nx_int32(this, j) result(nx)
  use, intrinsic :: iso_c_binding, only: c_long, c_int
  use atlas_grid_Structured_c_binding
  integer(c_long) :: nx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int), intent(in) :: j
  nx = atlas__grid__Structured__nx(this%c_ptr(), c_idx(j) )
end function

function Structured__nx_int64(this, j) result(nx)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  integer(c_long) :: nx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long), intent(in) :: j
  nx = atlas__grid__Structured__nx(this%c_ptr(), c_idx(j) )
end function

function Structured__reduced(this) result(reduced)
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in) :: this
  logical :: reduced
  if( atlas__grid__Structured__reduced(this%c_ptr()) == 1 ) then
    reduced = .true.
  else
    reduced = .false.
  endif
end function

function Structured__nx_array(this) result(nx)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_long, c_ptr, c_size_t, c_f_pointer
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long), pointer                :: nx(:)
  type   (c_ptr)                          :: nx_c_ptr
  integer(c_size_t)                       :: nx_size
  call atlas__grid__Structured__nx_array(this%c_ptr(), nx_c_ptr, nx_size)
  call c_f_pointer(nx_c_ptr , nx , (/nx_size/))
end function

function Structured__nxmax(this) result(nxmax)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  integer(c_long)                          :: nxmax
  nxmax = atlas__grid__Structured__nxmax(this%c_ptr())
end function

function Structured__nxmin(this) result(nxmin)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  integer(c_long)                          :: nxmin
  nxmin = atlas__grid__Structured__nxmin(this%c_ptr())
end function

function Structured__y(this, j) result(y)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: y
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long),             intent(in) :: j
  y = atlas__grid__Structured__y(this%c_ptr(), c_idx(j) )
end function

function Structured__y_32(this, j) result(y)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  real(c_double) :: y
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int),              intent(in) :: j
  y = atlas__grid__Structured__y(this%c_ptr(), c_idx(j) )
end function

function Structured__y_64(this, j) result(y)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: y
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long),             intent(in) :: j
  y = atlas__grid__Structured__y(this%c_ptr(), c_idx(j) )
end function

function Structured__y_array(this) result(y)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_double, c_ptr, c_size_t, c_f_pointer
  class(atlas_StructuredGrid), intent(in) :: this
  real   (c_double)       , pointer    :: y(:)
  type   (c_ptr)                       :: y_c_ptr
  integer(c_size_t)                    :: y_size
  call atlas__grid__Structured__y_array(this%c_ptr(), &
      & y_c_ptr, y_size)
  call c_f_pointer (y_c_ptr, y, (/y_size/))
end function

function Structured__x_32(this, i,j) result(x)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  real(c_double) :: x
  integer(c_int) :: i,j
  x = atlas__grid__Structured__x(this%c_ptr(), c_idx(i), c_idx(j))
end function

function Structured__x_64(this, i,j) result(x)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  real(c_double) :: x
  integer(c_long) :: i,j
  x = atlas__grid__Structured__x(this%c_ptr(), c_idx(i), c_idx(j))
end function

function Structured__xy_32(this, i,j) result(xy)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  real(c_double) :: xy(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int) , intent(in) :: i,j
  call atlas__grid__Structured__xy(this%c_ptr(), c_idx(i), c_idx(j), xy)
end function

function Structured__xy_64(this, i,j) result(xy)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: xy(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) , intent(in) :: i,j
  call atlas__grid__Structured__xy(this%c_ptr(), c_idx(i), c_idx(j), xy)
end function

function Structured__lonlat_32(this, i,j) result(lonlat)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  real(c_double) :: lonlat(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int) , intent(in) :: i,j
  call atlas__grid__Structured__lonlat(this%c_ptr(), c_idx(i), c_idx(j), lonlat)
end function

function Structured__lonlat_64(this, i,j) result(lonlat)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: lonlat(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) , intent(in) :: i,j
  call atlas__grid__Structured__lonlat(this%c_ptr(), c_idx(i), c_idx(j), lonlat)
end function

subroutine atlas_Grid__delete(this)
  use atlas_grid_Structured_c_binding
  class(atlas_Grid), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__grid__Structured__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine



subroutine atlas_Grid__copy(this,obj_in)
  class(atlas_Grid), intent(inout) :: this
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_Grid_module
