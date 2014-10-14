! (C) Copyright 2013-2014 ECMWF.


! -----------------------------------------------------------------------------
! ReducedGG routines

function new_reduced_gaussian_grid(identifier) result(rgg)
  type(ReducedGG_type) :: rgg
  character(len=*) :: identifier
  rgg%private%object = atlas__new_reduced_gaussian_grid(c_str(identifier))
end function new_reduced_gaussian_grid

function new_regular_gaussian_grid(nlon,nlat) result(rgg)
  type(ReducedGG_type) :: rgg
  integer, intent(in) :: nlon, nlat
  rgg%private%object = atlas__new_regular_gaussian_grid(nlon,nlat)
end function new_regular_gaussian_grid

function new_custom_reduced_gaussian_grid(nlon) result(rgg)
  type(ReducedGG_type) :: rgg
  integer, intent(in) :: nlon(:)
  rgg%private%object = atlas__new_custom_reduced_gaussian_grid(nlon,size(nlon))
end function new_custom_reduced_gaussian_grid

function new_regular_latlon_grid(nlon,nlat) result(rgg)
  type(ReducedGG_type) :: rgg
  integer, intent(in) :: nlon, nlat
  rgg%private%object = atlas__new_regular_latlon_grid(nlon,nlat)
end function new_regular_latlon_grid

function ReducedGG__ngptot(this) result(ngptot)
  class(ReducedGG_type), intent(in) :: this
  integer :: ngptot
  ngptot = atlas__RGG__ngptot(this%private%object)
end function ReducedGG__ngptot

function ReducedGG__nlat(this) result(nlat)
  class(ReducedGG_type), intent(in) :: this
  integer :: nlat
  nlat = atlas__RGG__nlat(this%private%object)
end function ReducedGG__nlat

function ReducedGG__nlon(this) result(nlon)
  class(ReducedGG_type), intent(in) :: this
  integer, pointer :: nlon(:)
  type(c_ptr) :: nlon_c_ptr
  integer(c_int) :: nlon_size
  call atlas__RGG__nlon(this%private%object, nlon_c_ptr, nlon_size)
  call C_F_POINTER ( nlon_c_ptr , nlon , (/nlon_size/) )
end function ReducedGG__nlon

function ReducedGG__lats(this) result(lats)
  class(ReducedGG_type), intent(in) :: this
  real(c_double), pointer :: lats(:)
  type(c_ptr) :: lats_c_ptr
  integer(c_int) :: lats_size
  call atlas__RGG__lats(this%private%object, lats_c_ptr, lats_size)
  call C_F_POINTER (  lats_c_ptr , lats , (/lats_size/) )
end function ReducedGG__lats

