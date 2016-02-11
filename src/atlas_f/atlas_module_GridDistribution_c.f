! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! GridDistribution routines

function atlas_GridDistribution__ctor( part, part0 ) result(griddistribution)
  type(atlas_GridDistribution) :: griddistribution
  integer, intent(in) :: part(:)
  integer, intent(in), optional :: part0
  integer:: npts, opt_part0
  opt_part0 = 0
  if( present(part0) ) opt_part0 = part0
  npts = size(part)
  call griddistribution%reset_c_ptr( atlas__GridDistribution__new(npts, part, opt_part0) );
end function atlas_GridDistribution__ctor


subroutine atlas_GridDistribution__final( this )
  class(atlas_GridDistribution), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__GridDistribution__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_GridDistribution__delete( this )
  class(atlas_GridDistribution), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__GridDistribution__delete(this%c_ptr());
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_GridDistribution__copy(this,obj_in)
  class(atlas_GridDistribution), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine


