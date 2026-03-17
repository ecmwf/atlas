! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_memory_resource
! This is a separate module which is used by pluto_module to implement the allocation and deallocation procedures which use the C memory resource interface.

use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, c_null_ptr, c_int32_t, c_int64_t, c_float, c_double, c_int, c_loc, &
                                      & c_f_pointer
implicit none
private

public :: pluto_memory_resource
public :: pluto_memory_pool_resource_reserve_int32
public :: pluto_memory_pool_resource_reserve_int64
public :: pluto_memory_pool_resource_reserve_real32
public :: pluto_memory_pool_resource_reserve_real64
public :: pluto_memory_pool_resource_release
public :: pluto_new_delete_resource
public :: pluto_null_memory_resource
public :: pluto_host_resource
public :: pluto_pinned_resource
public :: pluto_device_resource
public :: pluto_managed_resource
public :: pluto_mpi_resource
public :: pluto_host_pool_resource
public :: pluto_pinned_pool_resource
public :: pluto_device_pool_resource
public :: pluto_managed_pool_resource
public :: pluto_mpi_pool_resource
public :: pluto_has_registered_resource
public :: pluto_get_registered_resource
public :: pluto_register_resource
public :: pluto_unregister_resource
public :: pluto_set_label
public :: pluto_unset_label
public :: pluto_get_label
public :: pluto_register_memory_resource_adaptor

type :: pluto_memory_resource
    type(c_ptr) :: c_memory_resource = c_null_ptr
contains
    procedure :: allocate   => pluto_memory_resource_allocate
    procedure :: deallocate => pluto_memory_resource_deallocate
    procedure, private :: reserve_int32 => pluto_memory_pool_resource_reserve_int32
    procedure, private :: reserve_int64 => pluto_memory_pool_resource_reserve_int64
    procedure, private :: reserve_real32 => pluto_memory_pool_resource_reserve_real32
    procedure, private :: reserve_real64 => pluto_memory_pool_resource_reserve_real64
    generic, public :: reserve    => reserve_int32, reserve_int64, reserve_real32, reserve_real64
    procedure :: release    => pluto_memory_pool_resource_release
    procedure :: capacity   => pluto_memory_pool_resource_capacity
    procedure :: size       => pluto_memory_pool_resource_size
end type

contains

subroutine pluto_memory_resource_allocate(this, memory, bytes, alignment)
    class(pluto_memory_resource) :: this
    type(c_ptr), intent(out) :: memory
    integer(c_size_t), intent(in) :: bytes
    integer(c_size_t), intent(in), optional :: alignment
    interface
        function c_pluto_memory_resource_allocate(memory_resource, bytes, alignment) result(memory) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr) :: memory
            type(c_ptr), value :: memory_resource
            integer(c_size_t), value :: bytes
            integer(c_size_t), value :: alignment
        end function
    end interface

    if (present(alignment)) then
        memory = c_pluto_memory_resource_allocate(this%c_memory_resource, bytes, alignment)
    else
        memory = c_pluto_memory_resource_allocate(this%c_memory_resource, bytes, int(0,c_size_t))
    endif
end subroutine

subroutine pluto_memory_resource_deallocate(this, memory, bytes, alignment)
    class(pluto_memory_resource) :: this
    type(c_ptr), intent(inout) :: memory
    integer(c_size_t), intent(in) :: bytes
    integer(c_size_t), intent(in), optional :: alignment
    interface
        subroutine c_pluto_memory_resource_deallocate(memory_resource, memory, bytes, alignment) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value :: memory_resource
            type(c_ptr), value :: memory
            integer(c_size_t), value :: bytes
            integer(c_size_t), value :: alignment
        end subroutine
    end interface
    if (present(alignment)) then
        call c_pluto_memory_resource_deallocate(this%c_memory_resource, memory, bytes, alignment)
    else
        call c_pluto_memory_resource_deallocate(this%c_memory_resource, memory, bytes, int(0,c_size_t))
    endif
    memory = c_null_ptr
end subroutine


subroutine pluto_memory_pool_resource_release(this)
    class(pluto_memory_resource), intent(in) :: this
    interface
        subroutine c_pluto_memory_pool_resource_release(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr), value :: memory_resource
        end subroutine
    end interface
    call c_pluto_memory_pool_resource_release(this%c_memory_resource)
end subroutine

subroutine pluto_memory_pool_resource_reserve_size(this, bytes)
    class(pluto_memory_resource) :: this
    integer(c_size_t), intent(in) :: bytes
    interface
        subroutine c_pluto_memory_pool_resource_reserve(memory_resource, bytes) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value :: memory_resource
            integer(c_size_t), value :: bytes
        end subroutine
    end interface
    call c_pluto_memory_pool_resource_reserve(this%c_memory_resource, bytes)
end subroutine

subroutine pluto_memory_pool_resource_reserve_int32(this, bytes)
    class(pluto_memory_resource) :: this
    integer(c_int32_t), intent(in) :: bytes
    call pluto_memory_pool_resource_reserve_size(this, int(bytes,c_size_t))
end subroutine

subroutine pluto_memory_pool_resource_reserve_int64(this, bytes)
    class(pluto_memory_resource) :: this
    integer(c_int64_t), intent(in) :: bytes
    call pluto_memory_pool_resource_reserve_size(this, int(bytes,c_size_t))
end subroutine

subroutine pluto_memory_pool_resource_reserve_real32(this, bytes)
    class(pluto_memory_resource) :: this
    real(c_float), intent(in) :: bytes
    call pluto_memory_pool_resource_reserve_size(this, int(bytes,c_size_t))
end subroutine

subroutine pluto_memory_pool_resource_reserve_real64(this, bytes)
    class(pluto_memory_resource) :: this
    real(c_double), intent(in) :: bytes
    call pluto_memory_pool_resource_reserve_size(this, int(bytes,c_size_t))
end subroutine


function pluto_memory_pool_resource_size(this)
    integer(c_size_t) :: pluto_memory_pool_resource_size
    class(pluto_memory_resource), intent(in) :: this
    interface
        function c_pluto_memory_pool_resource_size(memory_resource) result(size) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            integer(c_size_t) :: size
            type(c_ptr), value :: memory_resource
        end function
    end interface
    pluto_memory_pool_resource_size = c_pluto_memory_pool_resource_size(this%c_memory_resource)
end function

function pluto_memory_pool_resource_capacity(this)
    integer(c_size_t) :: pluto_memory_pool_resource_capacity
    class(pluto_memory_resource), intent(in) :: this
    interface
        function c_pluto_memory_pool_resource_capacity(memory_resource) result(capacity) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            integer(c_size_t) :: capacity
            type(c_ptr), value :: memory_resource
        end function
    end interface
    pluto_memory_pool_resource_capacity = c_pluto_memory_pool_resource_capacity(this%c_memory_resource)
end function

function pluto_has_registered_resource(name)
    logical :: pluto_has_registered_resource
    character(len=*), target, intent(in) :: name
    integer(c_int) :: has_resource
    interface
        function c_pluto_has_registered_resource(name, name_size) result(has_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
            integer(c_int) :: has_resource
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
        end function
    end interface
    pluto_has_registered_resource = (c_pluto_has_registered_resource(c_loc(name), len_trim(name,kind=c_int)) /= 0)
end function

function pluto_get_registered_resource(name) result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    character(len=*), target, intent(in) :: name
    interface
        function c_pluto_get_registered_resource(name, name_size) result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
            type(c_ptr) :: memory_resource
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_get_registered_resource(c_loc(name), len_trim(name,kind=c_int))
end function

subroutine pluto_register_resource(name, memory_resource)
    character(len=*), target, intent(in) :: name
    type(pluto_memory_resource), intent(in) :: memory_resource
    interface
        subroutine c_pluto_register_resource(name, name_size, memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
            type(c_ptr), value, intent(in) :: memory_resource
        end subroutine
    end interface
    call c_pluto_register_resource(c_loc(name), len_trim(name,kind=c_int), memory_resource%c_memory_resource)
end subroutine

subroutine pluto_unregister_resource(name)
    character(len=*), target, intent(in) :: name
    interface
        subroutine c_pluto_unregister_resource(name, name_size) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
        end subroutine
    end interface
    call c_pluto_unregister_resource(c_loc(name), len_trim(name,kind=c_int))
end subroutine

function pluto_new_delete_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_new_delete_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_new_delete_resource()
end function

function pluto_null_memory_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_null_memory_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_null_memory_resource()
end function

function pluto_host_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_host_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_host_resource()
end function

function pluto_pinned_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_pinned_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_pinned_resource()
end function

function pluto_device_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_device_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_device_resource()
end function

function pluto_managed_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_managed_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_managed_resource()
end function

function pluto_mpi_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_mpi_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_mpi_resource()
end function

function pluto_pinned_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_pinned_pool_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_pinned_pool_resource()
end function

function pluto_host_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_host_pool_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_host_pool_resource()
end function

function pluto_device_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_device_pool_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_device_pool_resource()
end function

function pluto_managed_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_managed_pool_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_managed_pool_resource()
end function

function pluto_mpi_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    interface
        function c_pluto_mpi_pool_resource() result(memory_resource) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr
            type(c_ptr) :: memory_resource
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_mpi_pool_resource()
end function



subroutine pluto_set_label(label)
    implicit none
    character(len=*), target, intent(in) :: label
    interface
        subroutine c_pluto_set_label(label, label_size) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
            type(c_ptr), value, intent(in) :: label
            integer(c_int), value, intent(in) :: label_size
        end subroutine
    end interface
    call c_pluto_set_label(c_loc(label), len_trim(label))
end subroutine

subroutine pluto_unset_label()
    implicit none
    interface
        subroutine c_pluto_unset_label() bind(C)
        end subroutine
    end interface
    call c_pluto_unset_label()
end subroutine

function pluto_get_label() result(label)
    implicit none
    character(len=:), allocatable :: label
    type(c_ptr) :: label_cptr
    integer(c_int) :: label_size
    interface
        subroutine c_pluto_get_label(label, label_size) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
            type(c_ptr) :: label
            integer(c_int), intent(out) :: label_size
        end subroutine
    end interface
    call c_pluto_get_label(label_cptr, label_size)
    allocate(character(label_size) :: label)
    block
        character(len=1), pointer :: label_char_array(:)
        integer :: c
        call c_f_pointer(label_cptr, label_char_array, [label_size])
        do c=1,label_size
            label(c:c)=label_char_array(c)
        end do
    end block
end function

subroutine pluto_register_memory_resource_adaptor(name, allocate_fn, deallocate_fn)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_funptr, c_ptr, c_loc, c_funloc
    interface
        function allocate_signature(size, alignment) result(ptr) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr) :: ptr
            integer(c_size_t), value, intent(in) :: size
            integer(c_size_t), value, intent(in) :: alignment
        end function
        subroutine deallocate_signature(ptr, size, alignment) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value, intent(in) :: ptr
            integer(c_size_t), value, intent(in) :: size
            integer(c_size_t), value, intent(in) :: alignment
        end subroutine
        subroutine c_pluto_register_memory_resource_adaptor(name, name_size, allocate_fn, deallocate_fn) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_funptr
            type(c_ptr) :: memory_resource
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
            type(c_funptr), value :: allocate_fn
            type(c_funptr), value :: deallocate_fn
        end subroutine
    end interface
    character(len=*), target, intent(in) :: name
    procedure(allocate_signature)   :: allocate_fn
    procedure(deallocate_signature) :: deallocate_fn
    call c_pluto_register_memory_resource_adaptor(c_loc(name), len_trim(name,kind=c_int), &
        & c_funloc(allocate_fn), c_funloc(deallocate_fn))
end subroutine

end module
