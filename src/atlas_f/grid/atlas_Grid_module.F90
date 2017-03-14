
module atlas_Grid_module


use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_Grid
public :: atlas_StructuredGrid
public :: atlas_GaussianGrid
public :: atlas_ReducedGaussianGrid
public :: atlas_RegularGaussianGrid
public :: atlas_grid_RegularLonLat
public :: atlas_grid_ShiftedLonLat
public :: atlas_grid_ShiftedLon
public :: atlas_grid_ShiftedLat

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

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
  procedure :: npts => atlas_Grid__npts
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
  procedure :: nlat      => Structured__nlat
  procedure :: nlon      => Structured__nlon
  procedure :: pl        => Structured__pl
  procedure :: nlonmin   => Structured__nlonmin
  procedure :: nlonmax   => Structured__nlonmax
  procedure :: lat       => Structured__lat
  procedure :: latitudes => Structured__latitudes
  procedure :: lon       => Structured__lon
  procedure :: lonlat    => Structured__lonlat
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

TYPE, extends(atlas_StructuredGrid) :: atlas_grid_RegularLonLat

! Purpose :
! -------
!   *atlas_grid_RegularLonLat* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_RegularLonLat

interface atlas_grid_RegularLonLat
  module procedure atlas_grid_RegularLonLat__ctor_int32
  module procedure atlas_grid_RegularLonLat__ctor_int64
end interface


!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_grid_ShiftedLonLat

! Purpose :
! -------
!   *atlas_grid_ShiftedLonLat* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_ShiftedLonLat

interface atlas_grid_ShiftedLonLat
  module procedure atlas_grid_ShiftedLonLat__ctor_int32
  module procedure atlas_grid_ShiftedLonLat__ctor_int64
end interface
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_grid_ShiftedLon

! Purpose :
! -------
!   *atlas_grid_ShiftedLon* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_ShiftedLon

interface atlas_grid_ShiftedLon
  module procedure atlas_grid_ShiftedLon__ctor_int32
  module procedure atlas_grid_ShiftedLon__ctor_int64
end interface


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_grid_ShiftedLat

! Purpose :
! -------
!   *atlas_grid_ShiftedLat* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_ShiftedLat

interface atlas_grid_ShiftedLat
  module procedure atlas_grid_ShiftedLat__ctor_int32
  module procedure atlas_grid_ShiftedLat__ctor_int64
end interface

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!========================================================
contains
!========================================================

pure function c_idx(f_idx)
    use, intrinsic :: iso_c_binding, only : c_long
    integer(c_long) :: c_idx
    integer(c_long), intent(in) :: f_idx
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
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
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
  type(atlas_grid_RegularLonLat) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__RegularLonLat(int(nlon,c_long),int(nlat,c_long)) )
end function

function atlas_grid_RegularLonLat__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_RegularLonLat) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__RegularLonLat( nlon, nlat ) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLonLat__ctor_int32(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLonLat) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__ShiftedLonLat(int(nlon,c_long),int(nlat,c_long)) )
end function

function atlas_grid_ShiftedLonLat__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLonLat) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__ShiftedLonLat( nlon, nlat ) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLon__ctor_int32(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLon) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__ShiftedLon(int(nlon,c_size_t),int(nlat,c_size_t)) )
end function

function atlas_grid_ShiftedLon__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLon) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__ShiftedLon(int(nlon,c_size_t),int(nlat,c_size_t)) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLat__ctor_int32(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLat) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( &
    & atlas__grid__regular__ShiftedLat(int(nlon,c_size_t),int(nlat,c_size_t)))
end function

function atlas_grid_ShiftedLat__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLat) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( &
    & atlas__grid__regular__ShiftedLat(int(nlon,c_size_t),int(nlat,c_size_t)))
end function


! -----------------------------------------------------------------------------
! Structured members

function atlas_Grid__npts(this) result(npts)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_Grid), intent(in) :: this
  integer(c_long) :: npts
  npts = atlas__grid__Structured__npts(this%c_ptr())
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

function Structured__nlat(this) result(nlat)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) :: nlat
  nlat = atlas__grid__Structured__nlat(this%c_ptr())
end function


function Structured__nlon(this, jlat) result(nlon)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  integer(c_long) :: nlon
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long), intent(in) :: jlat
  nlon = atlas__grid__Structured__nlon(this%c_ptr(), c_idx(jlat) )
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

function Structured__pl(this) result(nlon)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_long, c_ptr, c_size_t, c_f_pointer
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long), pointer                 :: nlon(:)
  type   (c_ptr)                           :: nlon_c_ptr
  integer(c_size_t)                        :: nlon_size
  call atlas__grid__Structured__pl(this%c_ptr(), nlon_c_ptr, nlon_size)
  call c_f_pointer(nlon_c_ptr , nlon , (/nlon_size/))
end function

function Structured__nlonmax(this) result(nlonmax)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  integer(c_long)                           :: nlonmax
  nlonmax = atlas__grid__Structured__nlonmax(this%c_ptr())
end function

function Structured__nlonmin(this) result(nlonmin)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  integer(c_long)                           :: nlonmin
  nlonmin = atlas__grid__Structured__nlonmin(this%c_ptr())
end function

function Structured__lat(this, jlat) result(lat)
  use, intrinsic :: iso_c_binding, only: c_double, c_long, c_size_t
  use atlas_grid_Structured_c_binding
  real(c_double) :: lat
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long),              intent(in) :: jlat
  lat = atlas__grid__Structured__lat(this%c_ptr(), c_idx(jlat) )
end function

function Structured__latitudes(this) result(lat)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_double, c_ptr, c_size_t, c_f_pointer
  class(atlas_StructuredGrid), intent(in) :: this
  real   (c_double)       , pointer    :: lat(:)
  type   (c_ptr)                       :: lat_c_ptr
  integer(c_size_t)                    :: lat_size
  call atlas__grid__Structured__latitudes(this%c_ptr(), &
      & lat_c_ptr, lat_size)
  call c_f_pointer (lat_c_ptr, lat, (/lat_size/))
end function

function Structured__lon(this, jlat, jlon) result(lon)
  use, intrinsic :: iso_c_binding, only: c_double, c_long, c_size_t
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  real(c_double) :: lon
  integer(c_long) :: jlat
  integer(c_long) :: jlon
  lon = atlas__grid__Structured__lon(this%c_ptr(), c_idx(jlat), c_idx(jlon))
end function

function Structured__lonlat(this, jlat, jlon) result(lonlat)
  use, intrinsic :: iso_c_binding, only: c_double, c_long, c_size_t
  use atlas_grid_Structured_c_binding
  real(c_double) :: lonlat(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) , intent(in) :: jlat
  integer(c_long) , intent(in) :: jlon
  call atlas__grid__Structured__lonlat(this%c_ptr(), c_idx(jlat), c_idx(jlon), lonlat)
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
