module atlas_object_module
use iso_c_binding, only: c_ptr, c_null_ptr, c_associated
implicit none
private

private :: c_ptr
private :: c_null_ptr
private :: c_associated

!========================================================================
! Public interface

public atlas_object

!========================================================================

type, abstract :: atlas_object
  type(C_PTR), private :: cpp_object_ptr = C_NULL_PTR
contains
  procedure, public :: is_null => atlas_object__is_null
  procedure, public :: c_ptr   => atlas_object__c_ptr
  procedure, public :: reset_c_ptr => atlas_object__reset_c_ptr

  procedure, private :: equal => atlas_object__equal
  procedure, private :: not_equal => atlas_object__not_equal
  generic, public :: operator(==) => equal
  generic, public :: operator(/=) => not_equal

  procedure, public :: final => atlas_object__final
  procedure(atlas_object__final), deferred, public :: delete

end type

interface

  !int atlas__compare_cptr_equal( void* p1, void* p2 )
  function atlas__compare_cptr_equal(p1,p2) bind(c,name="atlas__compare_cptr_equal") result(equal)
    use iso_c_binding, only: c_ptr, c_int
    integer(c_int) :: equal
    type(c_ptr), value :: p1
    type(c_ptr), value :: p2
  end function

end interface


! =============================================================================
CONTAINS
! =============================================================================

function atlas_compare_equal(p1,p2) result(equal)
  use iso_c_binding, only: c_ptr
  logical :: equal
  type(c_ptr), intent(in) :: p1, p2
  if( atlas__compare_cptr_equal(p1,p2) == 1 ) then
    equal = .True.
  else
    equal = .False.
  endif
end function

function atlas_object__c_ptr(this) result(cptr)
  type(c_ptr) :: cptr
  class(atlas_object) :: this
  cptr = this%cpp_object_ptr
end function

subroutine atlas_object__reset_c_ptr(this,cptr)
  class(atlas_object) :: this
  type(c_ptr), optional :: cptr
  if( present(cptr) ) then
    this%cpp_object_ptr = cptr
  else
    this%cpp_object_ptr = c_null_ptr
  endif
end subroutine

function atlas_object__is_null(this) result(is_null)
  logical :: is_null
  class(atlas_object) :: this
  if( c_associated( this%cpp_object_ptr ) ) then
    is_null = .False.
  else
    is_null = .True.
  endif
end function

logical function atlas_object__equal(obj1,obj2)
  class(atlas_object), intent(in) :: obj1
  class(atlas_object), intent(in) :: obj2
  atlas_object__equal = atlas_compare_equal(obj1%c_ptr(),obj2%c_ptr())
end function

logical function atlas_object__not_equal(obj1,obj2)
  class(atlas_object), intent(in) :: obj1
  class(atlas_object), intent(in) :: obj2
  if( atlas_compare_equal(obj1%c_ptr(),obj2%c_ptr()) ) then
    atlas_object__not_equal = .False.
  else
    atlas_object__not_equal = .True.
  endif
end function


subroutine atlas_object__final(this)
  class(atlas_object), intent(inout) :: this
end subroutine

end module
