
program pluto_test_fortran_memory_resource

use pluto_module, only : pluto, pluto_memory_resource
use, intrinsic :: iso_c_binding, only: c_ptr, c_sizeof, c_float, c_f_pointer

implicit none

type(pluto_memory_resource) :: mr
call run_allocate(pluto%host%get_default_resource())
call run_allocate(pluto%new_delete_resource())
call run_allocate(pluto%managed_resource())
call run_allocate(pluto%pinned_resource())
call pluto%host%set_default_resource(pluto%managed_resource())
call run_allocate(pluto%host%get_default_resource())
call pluto%host%set_default_resource("pluto::pinned_resource")
call run_allocate(pluto%host%get_default_resource())

if (pluto%has_registered_resource("some_nonsense")) then
  write(0,*) "ERROR: resource 'some_nonsense' should not exist"
  ERROR STOP 1
endif

if (.NOT. pluto%has_registered_resource("pluto::host_resource")) then
  write(0,*) "ERROR: resource 'pluto::host_resource' should exist"
  ERROR STOP 1
endif

contains

subroutine run_allocate(memory_resource)
    implicit none
    integer, parameter :: wp = c_float
    type(pluto_memory_resource) :: memory_resource
    type(c_ptr) :: mem
    real(wp), pointer :: array1d(:)
    real(wp) :: real_value
    integer :: N = 10
    call memory_resource%allocate(mem, N*c_sizeof(real_value))
    call c_f_pointer(mem, array1d, [N])
    array1d(:) = 5
    call memory_resource%deallocate(mem, N*c_sizeof(real_value))
end subroutine

end program
