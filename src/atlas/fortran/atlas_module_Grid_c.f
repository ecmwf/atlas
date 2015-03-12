! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! ReducedGrid routines

function new_atlas_reduced_grid(identifier) result(grid)
  type(atlas_ReducedGrid) :: grid
  character(len=*) :: identifier
  grid%cpp_object_ptr = atlas__new_reduced_grid(c_str(identifier))
end function new_atlas_reduced_grid

function new_atlas_gaussian_grid(N) result(grid)
  type(atlas_ReducedGrid) :: grid
  integer, intent(in) :: N
  grid%cpp_object_ptr = atlas__new_gaussian_grid(N)
end function new_atlas_gaussian_grid

function new_atlas_reduced_gaussian_grid(nlon) result(grid)
  type(atlas_ReducedGrid) :: grid
  integer, intent(in) :: nlon(:)
  grid%cpp_object_ptr = atlas__new_reduced_gaussian_grid(nlon,size(nlon))
end function new_atlas_reduced_gaussian_grid

function new_atlas_lonlat_grid(nlon,nlat) result(grid)
  type(atlas_ReducedGrid) :: grid
  integer, intent(in) :: nlon, nlat
  grid%cpp_object_ptr = atlas__new_lonlat_grid(nlon,nlat)
end function new_atlas_lonlat_grid

function ReducedGrid__npts(this) result(npts)
  class(atlas_ReducedGrid), intent(in) :: this
  integer :: npts
  npts = atlas__ReducedGrid__npts(this%cpp_object_ptr)
end function ReducedGrid__npts

function ReducedGrid__nlat(this) result(nlat)
  class(atlas_ReducedGrid), intent(in) :: this
  integer :: nlat
  nlat = atlas__ReducedGrid__nlat(this%cpp_object_ptr)
end function ReducedGrid__nlat

function ReducedGrid__nlon(this) result(nlon)
  class(atlas_ReducedGrid), intent(in) :: this
  integer, pointer :: nlon(:)
  type(c_ptr) :: nlon_c_ptr
  integer(c_int) :: nlon_size
  call atlas__ReducedGrid__nlon(this%cpp_object_ptr, nlon_c_ptr, nlon_size)
  call C_F_POINTER ( nlon_c_ptr , nlon , (/nlon_size/) )
end function ReducedGrid__nlon

function ReducedGrid__latitudes(this) result(lat)
  class(atlas_ReducedGrid), intent(in) :: this
  real(c_double), pointer :: lat(:)
  type(c_ptr) :: lat_c_ptr
  integer(c_int) :: lat_size
  call atlas__ReducedGrid__latitudes(this%cpp_object_ptr, lat_c_ptr, lat_size)
  call C_F_POINTER (  lat_c_ptr , lat , (/lat_size/) )
end function ReducedGrid__latitudes

