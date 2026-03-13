! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------------------------------------------------
module my_allocator_mod
! ------------------------------------------------------------------------------------------------------------------------

use pluto_module, only : pluto_allocator, pluto_memory_resource
implicit none

type(pluto_memory_resource), save :: my_resource
    !! my_resource is a memory_resource which can be setup to any concrete implementation of a memory resource that is compatible with the pluto memory resource interface.
    !! It needs to be set up using my_allocator_init, which will also register it in pluto with the "my" string.
    !! In this example, we will set it to a mpi_pool memory pool resource, which is a predefined resource available in pluto that manages a memory pool with underlying
    !! (de)allocation using MPI_Alloc_mem / MPI_Free_mem.

type(pluto_allocator),       save :: my_allocator
    !! my_allocator is a pluto_allocator that uses my_resource as its memory resource.
    !! It needs to be set up using my_allocator_init.
    !! This allocator can then be used as a convenient API to allocate and deallocate memory using the my_resource, see examples below.

contains

subroutine my_allocator_init(bytes)
    !! Initialize my_resource and my_allocator module variables

    use pluto_module, only : pluto, pluto_memory_resource
    use iso_fortran_env, only : output_unit
    implicit none
    integer, intent(in) :: bytes ! the amount of memory to reserve initially in the pool resource, in bytes

    real(8), parameter :: GB = 1024**3 ! 1 GB in bytes

    write(0,'(A,I0,A)') "+ my_allocator_init(bytes=",bytes,")"

    call pluto%trace%enable() ! Just for this example to show the trace output of the resource,
                              ! can also be enabled by setting the environment variable PLUTO_TRACE=1

    ! Set my_resource to a memory pool that manages pinned host memory.
    ! This is a predefined resource that is available in pluto, and can be accessed directly as a member of the pluto object
    ! This resource will manage a pool of pinned host memory, which is useful for efficient data transfer between host and device.
    ! The resource will grow the pool as needed when allocations are made, but we can also reserve a certain amount of memory upfront
    ! to ensure that it is available when needed and to potentially improve performance by reducing fragmentation.
    my_resource = pluto%mpi_pool_resource()

    ! This can also be accessed by name, since it is registered with the pluto resource manager under the name "mpi_pool".
    !     my_resource = pluto%get_registered_resource("mpi_pool")

    ! Alternative could be a pinned memory pool resource:
    !     my_resource = pluto%pinned_pool_resource()
    !     my_resource = pluto%get_registered_resource("pinned_pool")

    ! NB: another implementation of a resource could easily be conceived, and registered by name, without having to change the rest of the code that uses it.
    ! For example using custom_resource_mod below:
    !    call register_custom_resource()
    !    my_resource = pluto%get_registered_resource("custom_resource") ! This is the resource that we registered with our custom allocate/deallocate functions.

    call my_resource%reserve(bytes) ! reserve memory if the resource supports it, otherwise this is a no-op

    ! Register my_resource by name in the pluto resource manager so that it can be accessed from other translation units
    ! without having to pass it around.
    call pluto%register_resource("MY", my_resource)

    ! Set my_allocator that uses my_resource as its memory resource.
    ! This allocator can then be used as a convenient API to allocate and deallocate memory using the my_resource.
    my_allocator = pluto%make_allocator(my_resource)
end subroutine my_allocator_init

subroutine my_allocator_finalize()
    use pluto_module, only : pluto
    write(0,'(A)') "+ my_allocator_finalize()"

    call my_resource%release()
    call pluto%unregister_resource("MY")

    ! Show memory report, listing statistics across all pluto-tracked memory resources
    ! This currently does not include "custom_resource", but would include the "MY" resource listed as "mpi" and "mpi_pool"
    write(0,'(A)') "Pluto Memory report:"
    write(0,'(A)') pluto%memory%report()
end subroutine my_allocator_finalize

subroutine my_allocator_print()
    use iso_fortran_env, only : output_unit
    implicit none
    real(8), parameter :: GB = 1024**3 ! 1 GB in bytes
    write(output_unit, '(A, F6.2, A, A, F6.2, A)') &
        "my_allocator: capacity:", my_resource%capacity()/GB, "GB", &
                  ",   allocated:", my_resource%size()/GB, "GB"
end subroutine my_allocator_print

end module my_allocator_mod
! ------------------------------------------------------------------------------------------------------------------------


! ------------------------------------------------------------------------------------------------------------------------
program main
! ------------------------------------------------------------------------------------------------------------------------

implicit none
call init_mpi()

! Initialise the my_allocator in a separate block to show that it can be done independently of the rest of the code,
! and that the resource is registered globally.
block
    use my_allocator_mod, only : my_allocator_init
    call my_allocator_init(1024**3) ! 1 GB
end block

! Example using my_allocator_mod, encapsulating pluto completely. This would be the recommended approach for the IFS.
block
    use my_allocator_mod, only : my_allocator
    real(8), pointer :: array(:,:)

    ! Allocation using shape, anonymous array
    call my_allocator%allocate(array, shape=[10, 8])
    call my_allocator%deallocate(array)

    ! Allocation using bounds, anonymous array
    call my_allocator%allocate(array, lbounds=[0,0], ubounds=[10, 8])
    call my_allocator%deallocate(array)

    ! label array for tracing
    call my_allocator%allocate("my_array", array, lbounds=[0,0], ubounds=[10, 8])
    call my_allocator%deallocate("my_array", array)
end block



! ------------------------------------------------------------------------------------------------------------------------
! Further examples that also use the pluto_module. This indicates that we will be able to use the same allocators
! in other contexts independent of MY, and that the resource is registered globally so that it can be used in other
! translation units without having to pass it around.

! Example using my_resource using pluto%allocate / pluto%deallocate
block
    use pluto_module, only : pluto
    use my_allocator_mod, only : my_resource
    real(8), pointer :: array(:,:)

    ! Allocation using shape, anonymous array
    call pluto%allocate(array, shape=[10, 8], resource=my_resource)
    call pluto%deallocate(array, resource=my_resource) ! IMPORTANT, must match the resource used for allocation

    ! Allocation using bounds, anonymous array
    call pluto%allocate(array, lbounds=[0,0], ubounds=[10, 8], resource=my_resource)
    call pluto%deallocate(array, resource=my_resource) ! IMPORTANT, must match the resource used for allocation

    ! label array for tracing
    call pluto%allocate("my_array_1", array, lbounds=[0,0], ubounds=[10, 8], resource=my_resource)
    call pluto%deallocate("my_array_1", array, resource=my_resource)

    ! Using label externally
    call pluto%set_label("my_array_2");
    call pluto%allocate(array, lbounds=[0,0], ubounds=[10, 8], resource=my_resource);
    call pluto%unset_label()
    !...
    call pluto%set_label("my_array_2");
    call pluto%deallocate(array, resource=my_resource);
    call pluto%unset_label()
end block

! Example independent of my_allocator_mod, using the resource name directly.
! This also shows that the resource is registered globally and can be used in other translation units.
block
    use pluto_module, only : pluto, pluto_allocator
    type(pluto_allocator) :: allocator
    real(8), pointer :: array(:,:)
    allocator = pluto%make_allocator("MY") ! This has been registered by name, now equivalent to using my_allocator

    ! anonymous array
    call allocator%allocate(array, lbounds=[0,0], ubounds=[10, 8])
    call allocator%deallocate(array)

    ! label array for tracing
    call allocator%allocate("my_array", array, lbounds=[0,0], ubounds=[10, 8])
    call allocator%deallocate("my_array", array)
end block

! Example independent of my_allocator_mod, using the resource name directly.
block
    use pluto_module, only : pluto
    real(8), pointer :: array(:,:)

    ! anonymous array
    call pluto%allocate(array, lbounds=[0,0], ubounds=[10, 8], resource="MY")
    call pluto%deallocate(array, resource="MY") ! IMPORTANT, must match the resource used for allocation

    ! label array for tracing
    call pluto%allocate("my_array", array, lbounds=[0,0], ubounds=[10, 8], resource="MY")
    call pluto%deallocate("my_array", array, resource="MY") ! IMPORTANT, must match the resource used for allocation

end block

! Example independent of my_allocator_mod, modifying the default host allocator in a scope.
! This is a convenient way to use the resource without having to pass it explicitly to every allocate/deallocate call.
block
    use pluto_module, only : pluto
    call pluto%scope%push()
    call pluto%host%set_default_resource("MY") ! set the default host resource for this scope to "MY"
    block
        real(8), pointer :: array(:,:)

        ! anonymous array
        call pluto%host%allocate(array, lbounds=[0,0], ubounds=[10, 8]) ! will use the default resource for the host, which is now "MY"
        call pluto%host%deallocate(array) ! will use the default resource for the host, which is now "MY"

        ! label array for tracing
        call pluto%host%allocate("my_array_1", array, lbounds=[0,0], ubounds=[10, 8]) ! will use the default resource for the host, which is now "MY"
        call pluto%host%deallocate("my_array_1", array) ! will use the default resource for the host, which is now "MY"

        ! Using shape argument instad of bounds
        call pluto%host%allocate("my_array_2", array, shape=[10, 8]) ! will use the default resource for the host, which is now my_resource
        call pluto%host%deallocate("my_array_2", array) ! will use the default resource for the host, which is now my_resource
    end block
    call pluto%scope%pop() ! restore the previous default resource for the host
end block

block
    use my_allocator_mod, only : my_allocator_finalize
    call my_allocator_finalize()
end block

call finalize_mpi()

contains

subroutine init_mpi()
    use pluto_module, only : pluto
    call pluto%mpi%init() ! calls MPI_Init, but only if mpirun is detected through typical environment variables
end subroutine

subroutine finalize_mpi()
    use pluto_module, only : pluto
    call pluto%mpi%finalize() ! calls MPI_Finalize, but only if mpirun is detected through typical environment variables
end subroutine


end program




!-----------------------------------------------------------------------------------------------------------------------------------------------
! Just to show how to register a custom memory resource using your own allocate and deallocate functions

module custom_resource_mod
public
contains

! This is an example of how a custom memory resource could be implemented and registered with pluto.
! It illustrates the flexibility of the pluto memory resource system,
! allowing users to define their own allocation strategies and integrate them seamlessly with the rest of the pluto ecosystem.
subroutine register_custom_resource()
    use pluto_module, only : pluto
    call pluto%register_memory_resource_adaptor("custom_resource", custom_allocate, custom_deallocate)
contains
    function custom_allocate(bytes, alignment) result(ptr) bind(C)
        use iso_c_binding, only: c_size_t, c_ptr, c_char, c_loc, c_null_ptr
        integer(c_size_t), value, intent(in) :: bytes
        integer(c_size_t), value, intent(in) :: alignment
        type(c_ptr) :: ptr
        ptr = c_null_ptr

        ! Here you would implement the actual allocation logic.
        block
            character(kind=c_char,len=1), pointer :: memory(:)
            if (bytes > 0) then
                allocate(memory(bytes))
                ptr = c_loc(memory(1))
                write(0,'(A,I0,A,I0,A,I0)') "custom_resource   allocated bytes=", bytes, " ptr=",loc(memory(1))
            else
                write(0,'(A,I0,A,I0,A,I0)') "custom_resource   allocated bytes=", bytes, " ptr=0"
            end if
        end block
    end function
    subroutine custom_deallocate(ptr, bytes, alignment) bind(C)
        use iso_c_binding, only: c_size_t, c_ptr, c_char, c_f_pointer, c_null_ptr
        type(c_ptr), value, intent(in) :: ptr
        integer(c_size_t), value, intent(in) :: bytes
        integer(c_size_t), value, intent(in) :: alignment

        ! Here you would implement the actual deallocation logic.
        block
            character(kind=c_char,len=1), pointer :: memory(:)
            if (bytes > 0) then
                call c_f_pointer(ptr, memory, [bytes])
                write(0,'(A,I0,A,I0,A,I0)') "custom_resource deallocated bytes=", bytes, " ptr=",loc(memory(1))
                deallocate(memory)
            else
                write(0,'(A,I0,A,I0,A,I0)') "custom_resource deallocated bytes=", bytes, " ptr=0"
            end if
        end block
    end subroutine
end subroutine register_custom_resource

end module
