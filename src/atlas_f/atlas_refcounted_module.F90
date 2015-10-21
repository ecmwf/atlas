module atlas_refcounted_module
use atlas_c_interop, only: atlas_compare_equal
use atlas_object_module, only: atlas_object
implicit none
private

private :: atlas_object
private :: atlas_compare_equal

!========================================================================
! Public interface

public atlas_RefCounted
public atlas_RefCounted_Fortran

!========================================================================

type, abstract, extends(atlas_object) :: atlas_RefCounted
contains
  procedure, public :: final => RefCounted__final
  procedure, private :: reset => RefCounted__reset
  generic, public :: assignment(=) => reset
  procedure(atlas_RefCounted__delete), deferred, public :: delete
  procedure, public :: owners => RefCounted__owners
  procedure, public :: attach => RefCounted__attach
  procedure, public :: detach => RefCounted__detach
  procedure, public :: return => atlas_RefCounted__return
endtype

interface
  subroutine atlas_RefCounted__delete(this)
     import atlas_RefCounted
     class(atlas_RefCounted), intent(inout):: this
  end subroutine
end interface


!========================================================================

type, abstract, extends(atlas_object) :: atlas_RefCounted_Fortran
  integer, private :: count = 0
contains
  procedure, public :: final => RefCounted_Fortran__final
  procedure, private :: reset => RefCounted_Fortran__reset
  generic, public :: assignment(=) => reset
  procedure(atlas_RefCounted_Fortran__delete), deferred, public :: delete
  procedure, public :: owners => RefCounted_Fortran__owners
  procedure, public :: attach => RefCounted_Fortran__attach
  procedure, public :: detach => RefCounted_Fortran__detach
  procedure, public :: return => atlas_RefCounted_Fortran__return
endtype

interface
  subroutine atlas_RefCounted_Fortran__delete(this)
     import atlas_RefCounted_Fortran
     class(atlas_RefCounted_Fortran), intent(inout):: this
  end subroutine
end interface

!========================================================================

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

!========================================================================
contains

subroutine RefCounted__final(this)
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
    if( .not. obj_out%is_null() ) call obj_out%final()
    call obj_out%reset_c_ptr( obj_in%c_ptr() )
    if( .not. obj_out%is_null() ) call obj_out%attach()
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

subroutine atlas_RefCounted__return(this)
  class(atlas_RefCounted), intent(inout) :: this
#ifndef FORTRAN_SUPPORTS_FINAL
  call this%detach()
#endif
end subroutine

!========================================================================

subroutine RefCounted_Fortran__final(this)
  class(atlas_RefCounted_Fortran), intent(inout) :: this
  if( .not. this%is_null() ) then
    if( this%owners() >  0 ) call this%detach()
    if( this%owners() == 0 ) call this%delete()
    call this%reset_c_ptr()
  endif
end subroutine

subroutine RefCounted_Fortran__reset(obj_out,obj_in)
  class(atlas_RefCounted_Fortran), intent(inout) :: obj_out
  class(atlas_RefCounted_Fortran), intent(in) :: obj_in
  if( obj_out /= obj_in ) then
    if( .not. obj_out%is_null() ) call obj_out%final()
    call obj_out%reset_c_ptr( obj_in%c_ptr() )
    obj_out%count = obj_in%count
    if( .not. obj_out%is_null() ) call obj_out%attach()
  endif
end subroutine

subroutine RefCounted_Fortran__attach(this)
  class(atlas_RefCounted_Fortran), intent(inout) :: this
  this%count = this%count + 1
end subroutine

subroutine RefCounted_Fortran__detach(this)
  class(atlas_RefCounted_Fortran), intent(inout) :: this
  this%count = this%count + 1
end subroutine

function RefCounted_Fortran__owners(this) result(owners)
  integer :: owners
  class(atlas_RefCounted_Fortran) :: this
  owners = this%count
end function

subroutine atlas_RefCounted_Fortran__return(this)
  class(atlas_RefCounted_Fortran), intent(inout) :: this
#ifndef FORTRAN_SUPPORTS_FINAL
  call this%detach()
#endif
end subroutine


end module
