program pluto_test_fortran_allocator

use pluto_module, only : pluto, pluto_allocator
use, intrinsic :: iso_c_binding, only : c_float

implicit none

type(pluto_allocator) :: allocator(2)

call run_allocate(pluto%host%make_allocator(), "program -- default")

call pluto%scope%push()
    call pluto%host%set_default_resource("pluto::pinned_resource")

    allocator(1) = pluto%host%make_allocator()
    call run_allocate(pluto%host%make_allocator(), "scope 1 -- pinned")

    call pluto%scope%push()
        call pluto%host%set_default_resource("pluto::managed_resource")

        allocator(2) = pluto%host%make_allocator()
        call run_allocate(pluto%host%make_allocator(), "scope 2 -- managed")

    call pluto%scope%pop()

    call run_allocate(pluto%host%make_allocator(), "scope 1 -- pinned")

call pluto%scope%pop()

call run_allocate(pluto%host%make_allocator(), "program -- default")

call run_allocate(allocator(1), "allocator(1) -- pinned")
call run_allocate(allocator(2), "allocator(2) -- managed")

contains

subroutine run_allocate(allocator, message)
    implicit none
    type(pluto_allocator) :: allocator
    character(len=*) :: message
    integer, parameter :: wp = c_float
    real(wp), pointer :: array4d(:,:,:,:)

    write(0,*) "run_allocate: ", message
    
    call allocator%allocate(array4d, [5,2,3,6])
    call allocator%deallocate(array4d)
end subroutine

end program
