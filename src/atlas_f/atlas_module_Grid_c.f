! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! ReducedGrid routines

function atlas_ReducedGrid__ctor_id(identifier) result(grid)
  type(atlas_ReducedGrid) :: grid
  character(len=*) :: identifier
  call grid%reset_c_ptr( atlas__new_reduced_grid(c_str(identifier)) )
end function atlas_ReducedGrid__ctor_id

function atlas_GaussianGrid__ctor(N) result(grid)
  type(atlas_GaussianGrid) :: grid
  integer, intent(in) :: N
  call grid%reset_c_ptr( atlas__new_gaussian_grid(N) )
end function atlas_GaussianGrid__ctor

function atlas_ReducedGaussianGrid__ctor(nlon) result(grid)
  type(atlas_ReducedGaussianGrid) :: grid
  integer, intent(in) :: nlon(:)
  call grid%reset_c_ptr( atlas__new_reduced_gaussian_grid(nlon,size(nlon)) )
end function atlas_ReducedGaussianGrid__ctor

function atlas_LonLatGrid__ctor(nlon,nlat) result(grid)
  type(atlas_LonLatGrid) :: grid
  integer, intent(in) :: nlon, nlat
  call grid%reset_c_ptr( atlas__new_lonlat_grid(nlon,nlat) )
end function atlas_LonLatGrid__ctor

function Grid__npts(this) result(npts)
  class(atlas_Grid), intent(in) :: this
  integer :: npts
  npts = atlas__ReducedGrid__npts(this%c_ptr())
end function Grid__npts

function ReducedGrid__N(this) result(N)
  class(atlas_ReducedGrid), intent(in) :: this
  integer :: N
  N = atlas__ReducedGrid__nlat(this%c_ptr())/2
end function ReducedGrid__N

function ReducedGrid__nlat(this) result(nlat)
  class(atlas_ReducedGrid), intent(in) :: this
  integer :: nlat
  nlat = atlas__ReducedGrid__nlat(this%c_ptr())
end function ReducedGrid__nlat

function ReducedGrid__nlon(this, jlat) result(nlon)
  class(atlas_ReducedGrid), intent(in) :: this
  integer                 , intent(in) :: jlat
  integer(8)                           :: nlon
  nlon = atlas__ReducedGrid__nlon(this%c_ptr(), jlat)
end function ReducedGrid__nlon

function ReducedGrid__nlon__all(this) result(nlon)
  class(atlas_ReducedGrid), intent(in) :: this
  integer, pointer                     :: nlon(:)
  type   (c_ptr)                       :: nlon_c_ptr
  integer(c_int)                       :: nlon_size
  call atlas__ReducedGrid__nlon__all(this%c_ptr(), nlon_c_ptr, nlon_size)
  call C_F_POINTER (nlon_c_ptr , nlon , (/nlon_size/))
end function ReducedGrid__nlon__all

function ReducedGrid__nlonmax(this) result(nlonmax)
  class(atlas_ReducedGrid), intent(in)  :: this
  integer                               :: nlonmax
  nlonmax = atlas__ReducedGrid__nlonmax(this%c_ptr())
end function ReducedGrid__nlonmax

function ReducedGrid__lat(this, jlat) result(lat)
  class(atlas_ReducedGrid), intent(in)  :: this
  real(c_double)                        :: lat
  integer                               :: jlat
  lat = atlas__ReducedGrid__lat(this%c_ptr(), jlat)
end function ReducedGrid__lat

function ReducedGrid__lat__all(this) result(lat)
  class(atlas_ReducedGrid), intent(in) :: this
  real   (c_double)       , pointer    :: lat(:)
  type   (c_ptr)                       :: lat_c_ptr
  integer(c_int)                       :: lat_size
  call atlas__ReducedGrid__lat__all(this%c_ptr(), lat_c_ptr, lat_size)
  call C_F_POINTER (lat_c_ptr, lat, (/lat_size/))
end function ReducedGrid__lat__all

function ReducedGrid__lon(this, jlat, jlon) result(lon)
  class(atlas_ReducedGrid), intent(in)  :: this
  real(c_double)                        :: lon
  integer                               :: jlat
  integer                               :: jlon
  lon = atlas__ReducedGrid__lon(this%c_ptr(), jlat, jlon)
end function ReducedGrid__lon

function ReducedGrid__lonlat(this, jlat, jlon) result(lonlat)
  class(atlas_ReducedGrid), intent(in) :: this
  integer                 , intent(in) :: jlat
  integer                 , intent(in) :: jlon
  real(c_double)          , pointer    :: lonlat(:)
  call atlas__ReducedGrid__lonlat(this%c_ptr(), jlat, jlon, lonlat)
end function ReducedGrid__lonlat

subroutine Grid__delete(this)
  class(atlas_Grid), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__ReducedGrid__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine Grid__delete

