! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! GridDistribution routines

function new_atlas_GridDistribution( part, part0 ) result(griddistribution)
  type(atlas_GridDistribution) :: griddistribution
  integer, intent(in) :: part(:)
  integer, intent(in), optional :: part0
  integer:: npts, opt_part0
  opt_part0 = 0
  if( present(part0) ) opt_part0 = part0
  npts = size(part)
  griddistribution%cpp_object_ptr = atlas__GridDistribution__new(npts, part, opt_part0);
end function new_atlas_GridDistribution


subroutine atlas_GridDistribution__delete( griddistribution )
  type(atlas_GridDistribution) :: griddistribution
  call atlas__GridDistribution__delete(griddistribution%cpp_object_ptr);
end subroutine




