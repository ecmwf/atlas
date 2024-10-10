program pluto_test_pluto_f

use pluto_module, only : pluto, pluto_memory_resource, pluto_allocator
use iso_c_binding

implicit none

integer, parameter :: wp = c_float

type(pluto_memory_resource) :: host_memory_resource
type(pluto_allocator) :: host_allocator
type(c_ptr) :: mem
real(wp), pointer :: array1d(:), array2d(:,:), array3d(:,:,:), array4d(:,:,:,:)
real(wp) :: real_value

call pluto%host%set_default_resource("pluto::managed_resource")
call pluto%host%get_default_resource(host_memory_resource)
call host_memory_resource%allocate(mem, 10*c_sizeof(real_value))

call c_f_pointer(mem, array1d, [10])
array1d(:) = 5
call host_memory_resource%deallocate(mem, 10*c_sizeof(real_value))

call pluto%host%get_default_allocator(host_allocator)

call host_allocator%allocate(array1d, [20])
call host_allocator%deallocate(array1d)

call host_allocator%allocate(array2d, [5,6])
call host_allocator%deallocate(array2d)

call host_allocator%allocate(array3d, [5,3,6])
call host_allocator%deallocate(array3d)

call host_allocator%allocate(array4d, [5,2,3,6])
call host_allocator%deallocate(array4d)

end program