! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! ReducedGrid routines

function new_reduced_gaussian_grid(identifier) result(grid)
  type(ReducedGrid_type) :: grid
  character(len=*) :: identifier
  grid%private%object = atlas__new_reduced_gaussian_grid(c_str(identifier))
end function new_reduced_gaussian_grid

function new_regular_gaussian_grid(nlon,nlat) result(grid)
  type(ReducedGrid_type) :: grid
  integer, intent(in) :: nlon, nlat
  grid%private%object = atlas__new_regular_gaussian_grid(nlon,nlat)
end function new_regular_gaussian_grid

function new_custom_reduced_gaussian_grid(nlon) result(grid)
  type(ReducedGrid_type) :: grid
  integer, intent(in) :: nlon(:)
  grid%private%object = atlas__new_custom_reduced_gaussian_grid(nlon,size(nlon))
end function new_custom_reduced_gaussian_grid

function new_regular_latlon_grid(nlon,nlat) result(grid)
  type(ReducedGrid_type) :: grid
  integer, intent(in) :: nlon, nlat
  grid%private%object = atlas__new_regular_latlon_grid(nlon,nlat)
end function new_regular_latlon_grid

function ReducedGrid__npts(this) result(npts)
  class(ReducedGrid_type), intent(in) :: this
  integer :: npts
  npts = atlas__ReducedGrid__npts(this%private%object)
end function ReducedGrid__npts

function ReducedGrid__nlat(this) result(nlat)
  class(ReducedGrid_type), intent(in) :: this
  integer :: nlat
  nlat = atlas__ReducedGrid__nlat(this%private%object)
end function ReducedGrid__nlat

function ReducedGrid__nlon(this) result(nlon)
  class(ReducedGrid_type), intent(in) :: this
  integer, pointer :: nlon(:)
  type(c_ptr) :: nlon_c_ptr
  integer(c_int) :: nlon_size
  call atlas__ReducedGrid__nlon(this%private%object, nlon_c_ptr, nlon_size)
  call C_F_POINTER ( nlon_c_ptr , nlon , (/nlon_size/) )
end function ReducedGrid__nlon

function ReducedGrid__latitudes(this) result(lat)
  class(ReducedGrid_type), intent(in) :: this
  real(c_double), pointer :: lat(:)
  type(c_ptr) :: lat_c_ptr
  integer(c_int) :: lat_size
  call atlas__ReducedGrid__latitudes(this%private%object, lat_c_ptr, lat_size)
  call C_F_POINTER (  lat_c_ptr , lat , (/lat_size/) )
end function ReducedGrid__latitudes

