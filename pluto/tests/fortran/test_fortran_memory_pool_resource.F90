program pluto_test_fortran_memory_pool_resource

use pluto_module, only : pluto, pluto_memory_resource, pluto_allocator
use iso_c_binding
implicit none

write(0,*) "reserve"
call reserve_mem()
write(0,*) "release"
call release_mem()

write(0,*) "allocate"
call run_allocate()
write(0,*) "release"
call release_mem()
write(0,*) "reserve"
call pluto%reserve(pluto%managed_pool_resource(), int(1024*1024*500,c_size_t))
write(0,*) "allocate"
call run_allocate()

! following should be a no-op, as new_delete_resource is not a memory_pool_resource
call pluto%reserve(pluto%new_delete_resource(), int(1024*1024*500,c_size_t))
call pluto%release(pluto%new_delete_resource())

write(0,*) "end"

contains

subroutine reserve_mem()
  type(pluto_memory_resource) :: memory_resource
  memory_resource = pluto%managed_pool_resource()
  call memory_resource%reserve(int(1024*1024*8,c_size_t))

  ! Or one-liner
  call pluto%reserve(pluto%managed_pool_resource(), int(1024*1024*257,c_size_t))
end subroutine

subroutine release_mem()
  type(pluto_memory_resource) :: memory_resource
  memory_resource = pluto%managed_pool_resource()
  call memory_resource%release()

  ! Or one-liner
  call pluto%release(pluto%managed_pool_resource())
end subroutine

subroutine run_allocate()
    implicit none
    integer :: j
    integer, parameter :: wp = c_float
    real(wp), pointer :: array2d(:,:)
    type(pluto_allocator) :: allocator
    allocator = pluto%make_allocator(pluto%managed_pool_resource())
    do j=1,100
        call allocator%allocate(array2d, shape=[1000*j, 1000])
        call allocator%deallocate(array2d)
    enddo
end subroutine

end program
