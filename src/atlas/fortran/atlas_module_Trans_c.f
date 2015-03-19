! (C) Copyright 2013-2015 ECMWF.


! -----------------------------------------------------------------------------
! Trans routines

#ifdef ATLAS_HAVE_TRANS
#define USE_ATLAS_TRANS_C_BINDING   use atlas_trans_c_binding
#else
#define USE_ATLAS_TRANS_C_BINDING
#define THROW_ERROR call atlas_throw_usererror("Cannot use atlas_Trans since atlas is compiled without ENABLE_TRANS=ON",atlas_code_location(__FILE__,__LINE__))
#endif

function new_atlas_Trans( grid ) result(trans)
  USE_ATLAS_TRANS_C_BINDING
  type(atlas_Trans) :: trans
  type(atlas_ReducedGrid) :: grid
#ifdef ATLAS_HAVE_TRANS
  trans%cpp_object_ptr = atlas__Trans__new( grid%cpp_object_ptr )
#else
  THROW_ERROR
#endif
end function new_atlas_Trans


subroutine atlas_Trans__delete( trans )
  USE_ATLAS_TRANS_C_BINDING
  type(atlas_Trans) :: trans
#ifdef ATLAS_HAVE_TRANS
  call atlas__Trans__delete(trans%cpp_object_ptr);
#else
  THROW_ERROR
#endif
end subroutine

function atlas_Trans__handle( this ) result(handle)
  USE_ATLAS_TRANS_C_BINDING
  integer :: handle
  class(atlas_Trans) :: this
#ifdef ATLAS_HAVE_TRANS
  handle = atlas__Trans__handle (this%cpp_object_ptr)
#else
  THROW_ERROR
#endif
end function



