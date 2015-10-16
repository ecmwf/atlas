module atlas_refcounted_module
use atlas_c_interop, only: atlas_compare_equal
use atlas_object_module, only: atlas_object
implicit none
private

private :: atlas_object
private :: atlas_compare_equal

!========================================================================
! Public interface

public atlas_refcounted

!========================================================================

type, abstract, extends(atlas_object) :: atlas_RefCounted
  integer, public :: id
contains
  procedure, public :: finalize => RefCounted__finalize
  procedure, private :: reset => RefCounted__reset
  generic, public :: assignment(=) => reset
  procedure(RefCounted__finalize), deferred, private :: delete
  procedure, public :: owners => RefCounted__owners
  procedure, public :: attach => RefCounted__attach
  procedure, public :: detach => RefCounted__detach
  procedure, public :: return => atlas_return
endtype

interface
  ! int atlas__Owned__owners(const Owned* This);
  function atlas__Owned__owners(This) bind(c,name="atlas__Owned__owners")
    use iso_c_binding, only: c_int, c_ptr
    integer(c_int) :: atlas__Owned__owners
    type(c_ptr), value :: This
  end function

  subroutine atlas__Owned__attach(This) bind(c,name="atlas__Owned__attach")
    use iso_c_binding, only: c_ptr
    type(c_ptr), value :: This
  end subroutine

  subroutine atlas__Owned__detach(This) bind(c,name="atlas__Owned__detach")
    use iso_c_binding, only: c_ptr
    type(c_ptr), value :: This
  end subroutine

end interface

contains

subroutine RefCounted__finalize(this)
  class(atlas_RefCounted), intent(inout) :: this
  if( .not. this%is_null() ) then
    if( this%owners() >  0 ) call this%detach()
    if( this%owners() == 0 ) call this%delete()
    call this%reset_c_ptr()
  endif
end subroutine

subroutine RefCounted__reset(obj_out,obj_in)
  class(atlas_RefCounted), intent(inout) :: obj_out
  class(atlas_RefCounted), intent(in) :: obj_in
  if( obj_out /= obj_in ) then
    if( .not. obj_out%is_null() ) call obj_out%finalize()
    call obj_out%reset_c_ptr( obj_in%c_ptr() )
    if( .not. obj_out%is_null() ) call obj_out%attach()
    obj_out%id = obj_in%id
  endif
end subroutine

subroutine RefCounted__attach(this)
  class(atlas_RefCounted), intent(inout) :: this
  call atlas__Owned__attach(this%c_ptr())
end subroutine

subroutine RefCounted__detach(this)
  class(atlas_RefCounted), intent(inout) :: this
  call atlas__Owned__detach(this%c_ptr())
end subroutine

function RefCounted__owners(this) result(owners)
  integer :: owners
  class(atlas_RefCounted) :: this
  owners = atlas__Owned__owners(this%c_ptr())
end function

subroutine atlas_return(this)
  class(atlas_RefCounted), intent(inout) :: this
#ifndef FORTRAN_SUPPORTS_FINAL
  call this%detach()
#endif
end subroutine


end module
