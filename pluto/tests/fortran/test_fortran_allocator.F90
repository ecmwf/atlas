program pluto_test_fortran_allocator

use pluto_module, only : pluto, pluto_allocator
use, intrinsic :: iso_c_binding, only : c_float

implicit none

call run_allocate( pluto%host%make_allocator() )
call pluto%host%set_default_resource("pluto::pinned_resource")
call run_allocate(pluto%host%make_allocator())
call run_allocate(pluto%make_allocator(pluto%managed_resource()))
call run_allocate(pluto%make_allocator("pluto::managed_resource"))
call run_allocate(pluto%make_allocator("pluto::managed_resource"))
call pluto%host%set_default_resource("pluto::pinned_resource")

contains

subroutine run_allocate(allocator)
    implicit none
    type(pluto_allocator) :: allocator
    integer, parameter :: wp = c_float
    real(wp), pointer :: array1d(:), array2d(:,:), array3d(:,:,:), array4d(:,:,:,:)

    call allocator%allocate(array1d, [20])
    call allocator%deallocate(array1d)
    
    call allocator%allocate(array2d, [5,6])
    call allocator%deallocate(array2d)
    
    call allocator%allocate(array3d, [5,3,6])
    call allocator%deallocate(array3d)
    
    call allocator%allocate(array4d, [5,2,3,6])
    call allocator%deallocate(array4d)
end subroutine

end program
