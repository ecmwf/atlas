
module atlas_Grid_module


use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_size_t, c_f_pointer
use atlas_c_interop, only: c_str
use atlas_refcounted_module, only: atlas_refcounted

implicit none

private :: c_ptr, c_int, c_double, c_size_t, c_f_pointer
private :: c_str
private :: atlas_refcounted

public :: atlas_Grid
public :: atlas_grid_Structured
public :: atlas_grid_ReducedGaussian
public :: atlas_grid_RegularGaussian
public :: atlas_grid_ShiftedLonLat

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_Grid

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

TYPE, extends(atlas_Grid) :: atlas_grid_Structured

! Purpose :
! -------
!   *atlas_grid_Structured* : Object Grid specifications for Reduced Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: N        => ReducedGrid__N
  procedure :: nlat     => ReducedGrid__nlat
  procedure :: nlon_idx => ReducedGrid__nlon
  procedure :: nlon_all => ReducedGrid__nlon__all
  generic   :: nlon     => nlon_idx, nlon_all
  procedure :: nlonmax  => ReducedGrid__nlonmax
  procedure :: lat_idx  => ReducedGrid__lat
  procedure :: lat_all  => ReducedGrid__lat__all
  generic   :: lat      => lat_idx, lat_all
  procedure :: lon      => ReducedGrid__lon
  procedure :: lonlat   => ReducedGrid__lonlat
END TYPE atlas_grid_Structured

!------------------------------------------------------------------------------

TYPE, extends(atlas_grid_Structured) :: atlas_grid_ReducedGaussian

! Purpose :
! -------
!   *atlas_grid_ReducedGaussian* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_ReducedGaussian

!------------------------------------------------------------------------------

TYPE, extends(atlas_grid_Structured) :: atlas_grid_RegularGaussian

! Purpose :
! -------
!   *atlas_grid_RegularGaussian* : Object Grid specifications for Regular Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_RegularGaussian

!------------------------------------------------------------------------------

TYPE, extends(atlas_grid_Structured) :: atlas_grid_ShiftedLonLat

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

!------------------------------------------------------------------------------

interface atlas_grid_Structured
  module procedure atlas_grid_Structured__ctor_id
  module procedure atlas_grid_Structured__ctor
end interface

interface atlas_grid_RegularGaussian
  module procedure atlas_grid_RegularGaussian__ctor
end interface

interface atlas_grid_ReducedGaussian
  module procedure atlas_grid_ReducedGaussian__ctor
end interface

interface atlas_grid_ShiftedLonLat
  module procedure atlas_grid_ShiftedLonLat__ctor
end interface

!------------------------------------------------------------------------------
!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! ReducedGrid routines

function atlas_grid_Structured__ctor_id(identifier) result(grid)
  use atlas_grid_global_Structured_c_binding
  type(atlas_grid_Structured) :: grid
  character(len=*), intent(in) :: identifier
  call grid%reset_c_ptr( atlas__new_reduced_grid(c_str(identifier)) )
end function atlas_grid_Structured__ctor_id

function atlas_grid_Structured__ctor(lats,nlon) result(grid)
  use atlas_grid_global_Structured_c_binding
  type(atlas_grid_Structured) :: grid
  real(c_double), intent(in) :: lats(:)
  integer, intent(in) :: nlon(:)
  call grid%reset_c_ptr( atlas__ReducedGrid__constructor(size(nlon),lats,nlon) )
end function atlas_grid_Structured__ctor

function atlas_grid_RegularGaussian__ctor(N) result(grid)
  use atlas_grid_global_Structured_c_binding
  type(atlas_grid_RegularGaussian) :: grid
  integer, intent(in) :: N
  call grid%reset_c_ptr( atlas__new_gaussian_grid(N) )
end function atlas_grid_RegularGaussian__ctor

function atlas_grid_ReducedGaussian__ctor(nlon) result(grid)
  use atlas_grid_global_Structured_c_binding
  type(atlas_grid_ReducedGaussian) :: grid
  integer, intent(in) :: nlon(:)
  call grid%reset_c_ptr( atlas__new_reduced_gaussian_grid(nlon,size(nlon)) )
end function atlas_grid_ReducedGaussian__ctor

function atlas_grid_ShiftedLonLat__ctor(nlon,nlat) result(grid)
  use atlas_grid_global_Structured_c_binding
  type(atlas_grid_ShiftedLonLat) :: grid
  integer, intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__new_lonlat_grid(nlon,nlat) )
end function atlas_grid_ShiftedLonLat__ctor

function atlas_Grid__npts(this) result(npts)
  use atlas_grid_global_Structured_c_binding
  class(atlas_Grid), intent(in) :: this
  integer :: npts
  npts = atlas__ReducedGrid__npts(this%c_ptr())
end function atlas_Grid__npts

function ReducedGrid__N(this) result(N)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer :: N
  N = atlas__ReducedGrid__nlat(this%c_ptr())/2
end function ReducedGrid__N

function ReducedGrid__nlat(this) result(nlat)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer :: nlat
  nlat = atlas__ReducedGrid__nlat(this%c_ptr())
end function ReducedGrid__nlat

function ReducedGrid__nlon(this, jlat) result(nlon)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer                 , intent(in) :: jlat
  integer(8)                           :: nlon
  nlon = atlas__ReducedGrid__nlon(this%c_ptr(), jlat-1)
end function ReducedGrid__nlon

function ReducedGrid__nlon__all(this) result(nlon)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer, pointer                     :: nlon(:)
  type   (c_ptr)                       :: nlon_c_ptr
  integer(c_int)                       :: nlon_size
  call atlas__ReducedGrid__nlon__all(this%c_ptr(), nlon_c_ptr, nlon_size)
  call C_F_POINTER (nlon_c_ptr , nlon , (/nlon_size/))
end function ReducedGrid__nlon__all

function ReducedGrid__nlonmax(this) result(nlonmax)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in)  :: this
  integer                               :: nlonmax
  nlonmax = atlas__ReducedGrid__nlonmax(this%c_ptr())
end function ReducedGrid__nlonmax

function ReducedGrid__lat(this, jlat) result(lat)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in)  :: this
  real(c_double)                        :: lat
  integer                               :: jlat
  lat = atlas__ReducedGrid__lat(this%c_ptr(), jlat-1)
end function ReducedGrid__lat

function ReducedGrid__lat__all(this) result(lat)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  real   (c_double)       , pointer    :: lat(:)
  type   (c_ptr)                       :: lat_c_ptr
  integer(c_int)                       :: lat_size
  call atlas__ReducedGrid__lat__all(this%c_ptr(), lat_c_ptr, lat_size)
  call C_F_POINTER (lat_c_ptr, lat, (/lat_size/))
end function ReducedGrid__lat__all

function ReducedGrid__lon(this, jlat, jlon) result(lon)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in)  :: this
  real(c_double)                        :: lon
  integer                               :: jlat
  integer                               :: jlon
  lon = atlas__ReducedGrid__lon(this%c_ptr(), jlat-1, jlon-1)
end function ReducedGrid__lon

function ReducedGrid__lonlat(this, jlat, jlon) result(lonlat)
  use atlas_grid_global_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer                 , intent(in) :: jlat
  integer                 , intent(in) :: jlon
  real(c_double)          , pointer    :: lonlat(:)
  call atlas__ReducedGrid__lonlat(this%c_ptr(), jlat-1, jlon-1, lonlat)
end function ReducedGrid__lonlat

subroutine atlas_Grid__delete(this)
  use atlas_grid_global_Structured_c_binding
  class(atlas_Grid), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__ReducedGrid__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_Grid__delete



subroutine atlas_Grid__copy(this,obj_in)
  class(atlas_Grid), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_Grid_module
