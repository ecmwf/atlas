
module atlas_Grid_module


use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_Grid
public :: atlas_grid_Structured
public :: atlas_grid_CustomStructured
public :: atlas_grid_ReducedGaussian
public :: atlas_grid_RegularGaussian
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
  procedure :: N         => Structured__N
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
END TYPE atlas_grid_Structured

interface atlas_grid_Structured
  module procedure atlas_grid_Structured__ctor_id
  module procedure atlas_grid_Structured__ctor_config
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_grid_Structured) :: atlas_grid_CustomStructured

! Purpose :
! -------
!   *atlas_grid_CustomStructured* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
END TYPE atlas_grid_CustomStructured

interface atlas_grid_CustomStructured
  module procedure atlas_grid_CustomStructured__ctor_int32
  module procedure atlas_grid_CustomStructured__ctor_int64
  module procedure atlas_grid_CustomStructured__ctor_lonmin_lonmax_int32
  module procedure atlas_grid_CustomStructured__ctor_lonmin_lonmax_int64
end interface

!------------------------------------------------------------------------------
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

interface atlas_grid_ReducedGaussian
  module procedure atlas_grid_ReducedGaussian__ctor_int32
  module procedure atlas_grid_ReducedGaussian__ctor_int64
end interface

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

interface atlas_grid_RegularGaussian
  module procedure atlas_grid_RegularGaussian__ctor_int32
  module procedure atlas_grid_RegularGaussian__ctor_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_grid_Structured) :: atlas_grid_RegularLonLat

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

interface atlas_grid_ShiftedLonLat
  module procedure atlas_grid_ShiftedLonLat__ctor_int32
  module procedure atlas_grid_ShiftedLonLat__ctor_int64
end interface
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

TYPE, extends(atlas_grid_Structured) :: atlas_grid_ShiftedLon

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

TYPE, extends(atlas_grid_Structured) :: atlas_grid_ShiftedLat

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
! -----------------------------------------------------------------------------
! Constructors

function atlas_grid_Structured__ctor_id(identifier) result(grid)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_grid_Structured) :: grid
  character(len=*), intent(in) :: identifier
  call grid%reset_c_ptr( atlas__grid__Structured(c_str(identifier)) )
end function

function atlas_grid_Structured__ctor_config(config) result(grid)
  use atlas_grid_Structured_c_binding
  use atlas_Config_module, only: atlas_Config
  type(atlas_grid_Structured) :: grid
  type(atlas_Config), intent(in) :: config
  call grid%reset_c_ptr( atlas__grid__Structured__config(config%c_ptr()) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_CustomStructured__ctor_int32(lats,nlon) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_CustomStructured) :: grid
  real(c_double), intent(in) :: lats(:)
  integer(c_int), intent(in) :: nlon(:)
  integer(c_size_t) :: nlat
  nlat = size(nlon)
  call grid%reset_c_ptr( atlas__grid__CustomStructured_int(nlat,lats,nlon) )
end function

function atlas_grid_CustomStructured__ctor_int64(lats,nlon) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_double, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_CustomStructured) :: grid
  real(c_double), intent(in) :: lats(:)
  integer(c_long), intent(in) :: nlon(:)
  integer(c_size_t) :: nlat
  nlat = size(nlon)
  call grid%reset_c_ptr( atlas__grid__CustomStructured_long(nlat,lats,nlon) )
end function

function atlas_grid_CustomStructured__ctor_lonmin_lonmax_int32(lats,nlon,lonmin,lonmax) result(grid)
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_CustomStructured) :: grid
  real(c_double), intent(in) :: lats(:)
  integer(c_int), intent(in) :: nlon(:)
  real(c_double), intent(in) :: lonmin(:)
  real(c_double), intent(in) :: lonmax(:)
  integer(c_size_t) :: nlat
  nlat = size(nlon)
  call grid%reset_c_ptr( atlas__grid__CustomStructured_lonmin_lonmax_int(nlat,lats,nlon,lonmin,lonmax) )
end function

function atlas_grid_CustomStructured__ctor_lonmin_lonmax_int64(lats,nlon,lonmin,lonmax) result(grid)
  use, intrinsic :: iso_c_binding, only: c_double, c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_CustomStructured) :: grid
  real(c_double),  intent(in) :: lats(:)
  integer(c_long), intent(in) :: nlon(:)
  real(c_double),  intent(in) :: lonmin(:)
  real(c_double),  intent(in) :: lonmax(:)
  integer(c_size_t) :: nlat
  nlat = size(nlon)
  call grid%reset_c_ptr( atlas__grid__CustomStructured_lonmin_lonmax_long(nlat,lats,nlon,lonmin,lonmax) )
end function



!-----------------------------------------------------------------------------

function atlas_grid_RegularGaussian__ctor_int32(N) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_RegularGaussian) :: grid
  integer(c_int), intent(in) :: N
  call grid%reset_c_ptr( atlas__grid__regular__RegularGaussian(int(N,c_size_t)) )
end function

function atlas_grid_RegularGaussian__ctor_int64(N) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_RegularGaussian) :: grid
  integer(c_long), intent(in) :: N
  call grid%reset_c_ptr( atlas__grid__regular__RegularGaussian(int(N,c_size_t)) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_ReducedGaussian__ctor_int32(N,nlon) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ReducedGaussian) :: grid
  integer(c_int), intent(in) :: N
  integer(c_int), intent(in)  :: nlon(:)
  call grid%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_int(int(N,c_size_t),nlon) )
end function

function atlas_grid_ReducedGaussian__ctor_int64(N,nlon) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ReducedGaussian) :: grid
  integer(c_long), intent(in) :: N
  integer(c_long), intent(in)  :: nlon(:)
  call grid%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_long(int(N,c_size_t),nlon) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_RegularLonLat__ctor_int32(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_RegularLonLat) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__RegularLonLat(int(nlon,c_size_t),int(nlat,c_size_t)) )
end function

function atlas_grid_RegularLonLat__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_RegularLonLat) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__RegularLonLat(int(nlon,c_size_t),int(nlat,c_size_t)) )
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLonLat__ctor_int32(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLonLat) :: grid
  integer(c_int), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__ShiftedLonLat(int(nlon,c_size_t),int(nlat,c_size_t)) )
end function

function atlas_grid_ShiftedLonLat__ctor_int64(nlon,nlat) result(grid)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t
  use atlas_grid_Structured_c_binding
  type(atlas_grid_ShiftedLonLat) :: grid
  integer(c_long), intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__grid__regular__ShiftedLonLat(int(nlon,c_size_t),int(nlat,c_size_t)) )
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

function Structured__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Structured__N(this%c_ptr())
end function

function Structured__nlat(this) result(nlat)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
  integer(c_long) :: nlat
  nlat = atlas__grid__Structured__nlat(this%c_ptr())
end function


function Structured__nlon(this, jlat) result(nlon)
  use, intrinsic :: iso_c_binding, only: c_long, c_size_t, c_int
  use atlas_grid_Structured_c_binding
  integer(c_long) :: nlon
  class(atlas_grid_Structured), intent(in) :: this
  integer(c_int), intent(in) :: jlat
  nlon = atlas__grid__Structured__nlon(this%c_ptr(), int(jlat-1,c_size_t) )
end function

function Structured__reduced(this) result(reduced)
  use atlas_grid_Structured_c_binding
  class(atlas_grid_Structured), intent(in) :: this
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
  class(atlas_grid_Structured), intent(in) :: this
  integer(c_long), pointer                 :: nlon(:)
  type   (c_ptr)                           :: nlon_c_ptr
  integer(c_size_t)                        :: nlon_size
  call atlas__grid__Structured__pl(this%c_ptr(), nlon_c_ptr, nlon_size)
  call c_f_pointer(nlon_c_ptr , nlon , (/nlon_size/))
end function

function Structured__nlonmax(this) result(nlonmax)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_grid_Structured), intent(in)  :: this
  integer(c_long)                           :: nlonmax
  nlonmax = atlas__grid__Structured__nlonmax(this%c_ptr())
end function

function Structured__nlonmin(this) result(nlonmin)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_grid_Structured), intent(in)  :: this
  integer(c_long)                           :: nlonmin
  nlonmin = atlas__grid__Structured__nlonmin(this%c_ptr())
end function

function Structured__lat(this, jlat) result(lat)
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_size_t
  use atlas_grid_Structured_c_binding
  real(c_double) :: lat
  class(atlas_grid_Structured), intent(in) :: this
  integer(c_int),               intent(in) :: jlat
  lat = atlas__grid__Structured__lat(this%c_ptr(), int(jlat-1,c_size_t))
end function

function Structured__latitudes(this) result(lat)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_double, c_ptr, c_size_t, c_f_pointer
  class(atlas_grid_Structured), intent(in) :: this
  real   (c_double)       , pointer    :: lat(:)
  type   (c_ptr)                       :: lat_c_ptr
  integer(c_size_t)                    :: lat_size
  call atlas__grid__Structured__latitudes(this%c_ptr(), &
      & lat_c_ptr, lat_size)
  call c_f_pointer (lat_c_ptr, lat, (/lat_size/))
end function

function Structured__lon(this, jlat, jlon) result(lon)
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_size_t
  use atlas_grid_Structured_c_binding
  class(atlas_grid_Structured), intent(in)  :: this
  real(c_double) :: lon
  integer(c_int) :: jlat
  integer(c_int) :: jlon
  lon = atlas__grid__Structured__lon(this%c_ptr(), &
      & int(jlat-1,c_size_t), int(jlon-1,c_size_t) )
end function

function Structured__lonlat(this, jlat, jlon) result(lonlat)
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_size_t
  use atlas_grid_Structured_c_binding
  real(c_double) :: lonlat(2)
  class(atlas_grid_Structured), intent(in) :: this
  integer(c_int) , intent(in) :: jlat
  integer(c_int) , intent(in) :: jlon
  call atlas__grid__Structured__lonlat(this%c_ptr(), &
      & int(jlat-1,c_size_t), int(jlon-1,c_size_t), lonlat)
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
