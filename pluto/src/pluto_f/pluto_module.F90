module pluto_module

use iso_c_binding, only: c_loc, c_ptr, c_int, c_size_t, c_null_ptr, c_double, c_float, c_int32_t, c_int64_t, c_f_pointer
implicit none
private

public :: pluto, pluto_memory_resource, pluto_allocator

interface
    function c_pluto_devices() result(devices) bind(c)
        use iso_c_binding, only: c_int
        integer(c_int) :: devices
    end function
    subroutine c_pluto_host_set_default_resource_name(name, name_size) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end subroutine
    subroutine c_pluto_host_set_default_resource_ptr(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr), value, intent(in) :: memory_resource
    end subroutine
    subroutine c_pluto_device_set_default_resource_name(name, name_size) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end subroutine
    subroutine c_pluto_device_set_default_resource_ptr(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr), value, intent(in) :: memory_resource
    end subroutine
    function c_pluto_host_get_default_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_device_get_default_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_has_registered_resource(name, name_size) result(has_resource) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        integer(c_int) :: has_resource
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end function
    function c_pluto_get_registered_resource(name, name_size) result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr) :: memory_resource
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end function
    subroutine c_pluto_register_resource(name, name_size, memory_resource) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
        type(c_ptr), value, intent(in) :: memory_resource
    end subroutine
    subroutine c_pluto_unregister_resource(name, name_size) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end subroutine
    function c_pluto_new_delete_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_null_memory_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_host_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_pinned_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_device_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_managed_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_host_pool_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_pinned_pool_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_device_pool_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_managed_pool_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_memory_resource_allocate(memory_resource, bytes, alignment) result(memory) bind(c)
        use iso_c_binding, only: c_ptr, c_size_t
        type(c_ptr) :: memory
        type(c_ptr), value :: memory_resource
        integer(c_size_t), value :: bytes
        integer(c_size_t), value :: alignment
    end function
    subroutine c_pluto_memory_resource_deallocate(memory_resource, memory, bytes, alignment) bind(c)
        use iso_c_binding, only: c_ptr, c_size_t
        type(c_ptr), value :: memory_resource
        type(c_ptr), value :: memory
        integer(c_size_t), value :: bytes
        integer(c_size_t), value :: alignment
    end subroutine
    function c_pluto_memory_pool_resource_size(memory_resource) result(size) bind(c)
        use iso_c_binding, only: c_ptr, c_size_t
        integer(c_size_t) :: size
        type(c_ptr), value :: memory_resource
    end function
    function c_pluto_memory_pool_resource_capacity(memory_resource) result(capacity) bind(c)
        use iso_c_binding, only: c_ptr, c_size_t
        integer(c_size_t) :: capacity
        type(c_ptr), value :: memory_resource
    end function
    subroutine c_pluto_memory_pool_resource_reserve(memory_resource, bytes) bind(c)
        use iso_c_binding, only: c_ptr, c_size_t
        type(c_ptr), value :: memory_resource
        integer(c_size_t), value :: bytes
    end subroutine
    subroutine c_pluto_memory_pool_resource_release(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr), value :: memory_resource
    end subroutine
    subroutine c_pluto_scope_push() bind(c)
    end subroutine
    subroutine c_pluto_scope_pop() bind(c)
    end subroutine
    subroutine c_pluto_trace_enable(enable) bind(c)
        use iso_c_binding, only: c_int
        integer(c_int), value :: enable
    end subroutine
    function c_pluto_trace_enabled() result(enabled) bind(c)
        use iso_c_binding, only: c_int
        integer(c_int) :: enabled
    end function
end interface

type :: pluto_trace_t
contains
    procedure, nopass :: enable  => pluto_trace_enable
    procedure, nopass :: enabled => pluto_trace_enabled
end type

type :: pluto_memory_resource
    type(c_ptr) :: c_memory_resource
contains
    procedure :: allocate   => pluto_memory_resource_allocate
    procedure :: deallocate => pluto_memory_resource_deallocate
    procedure :: reserve    => pluto_memory_pool_resource_reserve
    procedure :: release    => pluto_memory_pool_resource_release
    procedure :: capacity   => pluto_memory_pool_resource_capacity
    procedure :: size       => pluto_memory_pool_resource_size
end type

type pluto_allocator
    type(pluto_memory_resource) :: memory_resource
contains

    procedure :: pluto_allocator_allocate_int32_r1
    procedure :: pluto_allocator_allocate_int32_r2
    procedure :: pluto_allocator_allocate_int32_r3
    procedure :: pluto_allocator_allocate_int32_r4
    procedure :: pluto_allocator_allocate_int64_r1
    procedure :: pluto_allocator_allocate_int64_r2
    procedure :: pluto_allocator_allocate_int64_r3
    procedure :: pluto_allocator_allocate_int64_r4
    procedure :: pluto_allocator_allocate_real32_r1
    procedure :: pluto_allocator_allocate_real32_r2
    procedure :: pluto_allocator_allocate_real32_r3
    procedure :: pluto_allocator_allocate_real32_r4
    procedure :: pluto_allocator_allocate_real64_r1
    procedure :: pluto_allocator_allocate_real64_r2
    procedure :: pluto_allocator_allocate_real64_r3
    procedure :: pluto_allocator_allocate_real64_r4

    generic :: allocate => &
        & pluto_allocator_allocate_int32_r1, &
        & pluto_allocator_allocate_int32_r2, &
        & pluto_allocator_allocate_int32_r3, &
        & pluto_allocator_allocate_int32_r4, &
        & pluto_allocator_allocate_int64_r1, &
        & pluto_allocator_allocate_int64_r2, &
        & pluto_allocator_allocate_int64_r3, &
        & pluto_allocator_allocate_int64_r4, &
        & pluto_allocator_allocate_real32_r1, &
        & pluto_allocator_allocate_real32_r2, &
        & pluto_allocator_allocate_real32_r3, &
        & pluto_allocator_allocate_real32_r4, &
        & pluto_allocator_allocate_real64_r1, &
        & pluto_allocator_allocate_real64_r2, &
        & pluto_allocator_allocate_real64_r3, &
        & pluto_allocator_allocate_real64_r4

    procedure :: pluto_allocator_deallocate_int32_r1
    procedure :: pluto_allocator_deallocate_int32_r2
    procedure :: pluto_allocator_deallocate_int32_r3
    procedure :: pluto_allocator_deallocate_int32_r4
    procedure :: pluto_allocator_deallocate_int64_r1
    procedure :: pluto_allocator_deallocate_int64_r2
    procedure :: pluto_allocator_deallocate_int64_r3
    procedure :: pluto_allocator_deallocate_int64_r4
    procedure :: pluto_allocator_deallocate_real32_r1
    procedure :: pluto_allocator_deallocate_real32_r2
    procedure :: pluto_allocator_deallocate_real32_r3
    procedure :: pluto_allocator_deallocate_real32_r4
    procedure :: pluto_allocator_deallocate_real64_r1
    procedure :: pluto_allocator_deallocate_real64_r2
    procedure :: pluto_allocator_deallocate_real64_r3
    procedure :: pluto_allocator_deallocate_real64_r4

    generic :: deallocate => &
        & pluto_allocator_deallocate_int32_r1, &
        & pluto_allocator_deallocate_int32_r2, &
        & pluto_allocator_deallocate_int32_r3, &
        & pluto_allocator_deallocate_int32_r4, &
        & pluto_allocator_deallocate_int64_r1, &
        & pluto_allocator_deallocate_int64_r2, &
        & pluto_allocator_deallocate_int64_r3, &
        & pluto_allocator_deallocate_int64_r4, &
        & pluto_allocator_deallocate_real32_r1, &
        & pluto_allocator_deallocate_real32_r2, &
        & pluto_allocator_deallocate_real32_r3, &
        & pluto_allocator_deallocate_real32_r4, &
        & pluto_allocator_deallocate_real64_r1, &
        & pluto_allocator_deallocate_real64_r2, &
        & pluto_allocator_deallocate_real64_r3, &
        & pluto_allocator_deallocate_real64_r4
end type

type pluto_host_t
   integer :: dummy
contains
    procedure, nopass :: get_default_resource  => pluto_host_get_default_resource
    procedure, nopass :: make_allocator        => pluto_host_make_allocator
    procedure, private :: set_default_resource_name  => pluto_host_set_default_resource_name
    procedure, private :: set_default_resource_type  => pluto_host_set_default_resource_type
    generic, public :: set_default_resource => set_default_resource_type, set_default_resource_name
end type

type pluto_device_t
contains
    procedure, nopass :: get_default_resource  => pluto_device_get_default_resource
    procedure, nopass :: make_allocator        => pluto_device_make_allocator
    procedure, private :: set_default_resource_name  => pluto_device_set_default_resource_name
    procedure, private :: set_default_resource_type  => pluto_device_set_default_resource_type
    generic :: set_default_resource => set_default_resource_name, set_default_resource_type
end type

type pluto_scope_t
contains
    procedure, nopass :: push => pluto_scope_push
    procedure, nopass :: pop  => pluto_scope_pop
end type

type pluto_t
    type(pluto_host_t)   :: host
    type(pluto_device_t) :: device
    type(pluto_scope_t)  :: scope
    type(pluto_trace_t)  :: trace
contains
    procedure, nopass :: devices => pluto_devices
    procedure, nopass :: has_registered_resource => pluto_has_registered_resource
    procedure, nopass :: get_registered_resource => pluto_get_registered_resource
    procedure, nopass :: register_resource => pluto_register_resource
    procedure, nopass :: unregister_resource => pluto_unregister_resource
    procedure, nopass :: new_delete_resource => pluto_new_delete_resource
    procedure, nopass :: null_memory_resource => pluto_null_memory_resource
    procedure, nopass :: host_resource => pluto_host_resource
    procedure, nopass :: pinned_resource => pluto_pinned_resource
    procedure, nopass :: device_resource => pluto_device_resource
    procedure, nopass :: managed_resource => pluto_managed_resource
    procedure, nopass :: host_pool_resource => pluto_host_pool_resource
    procedure, nopass :: pinned_pool_resource => pluto_pinned_pool_resource
    procedure, nopass :: device_pool_resource => pluto_device_pool_resource
    procedure, nopass :: managed_pool_resource => pluto_managed_pool_resource
    procedure, private :: make_allocator_type => pluto_make_allocator_type
    procedure, private :: make_allocator_name => pluto_make_allocator_name
    generic :: make_allocator => make_allocator_type, make_allocator_name
    procedure, nopass :: reserve => pluto_memory_pool_resource_reserve
    procedure, nopass :: release => pluto_memory_pool_resource_release
end type

type(pluto_t) :: pluto

contains

subroutine pluto_host_set_default_resource_name(this, name)
    class(pluto_host_t) :: this
    character(len=*), target, intent(in) :: name
    call c_pluto_host_set_default_resource_name(c_loc(name), len(name,kind=c_int))
end subroutine

subroutine pluto_host_set_default_resource_type(this, memory_resource)
    class(pluto_host_t) :: this
    type(pluto_memory_resource), intent(in) :: memory_resource
    call c_pluto_host_set_default_resource_ptr(memory_resource%c_memory_resource)
end subroutine

subroutine pluto_device_set_default_resource_name(this, name)
    class(pluto_device_t) :: this
    character(len=*), target, intent(in) :: name
    call c_pluto_device_set_default_resource_name(c_loc(name), len(name,kind=c_int))
end subroutine

subroutine pluto_device_set_default_resource_type(this, memory_resource)
    class(pluto_device_t) :: this
    type(pluto_memory_resource), intent(in) :: memory_resource
    call c_pluto_device_set_default_resource_ptr(memory_resource%c_memory_resource)
end subroutine

subroutine pluto_scope_push()
    call c_pluto_scope_push()
end subroutine

subroutine pluto_scope_pop()
    call c_pluto_scope_pop()
end subroutine

subroutine pluto_trace_enable(enable)
    logical, optional :: enable
    logical :: do_enable
    do_enable = .true.
    if (present(enable)) then
        do_enable = enable
    endif
    if (do_enable) then
        call c_pluto_trace_enable(1_c_int)
    else
        call c_pluto_trace_enable(0_c_int)
    endif
end subroutine

function pluto_trace_enabled()
    logical :: pluto_trace_enabled
    integer(c_int) :: enabled
    enabled = c_pluto_trace_enabled()
    if (enabled == 1) then
        pluto_trace_enabled = .true.
    else
        pluto_trace_enabled = .false.
    endif
end function

function pluto_devices()
    integer(c_int) :: pluto_devices
    pluto_devices = c_pluto_devices()
end function

function pluto_has_registered_resource(name)
    logical :: pluto_has_registered_resource
    character(len=*), target, intent(in) :: name
    integer(c_int) :: has_resource
    has_resource = &
        & c_pluto_has_registered_resource(c_loc(name), len(name,kind=c_int))
end function

function pluto_get_registered_resource(name) result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    character(len=*), target, intent(in) :: name
    memory_resource%c_memory_resource = &
        & c_pluto_get_registered_resource(c_loc(name), len(name,kind=c_int))
end function

subroutine pluto_register_resource(name, memory_resource)
    character(len=*), target, intent(in) :: name
    type(pluto_memory_resource), intent(in) :: memory_resource
    call c_pluto_register_resource(c_loc(name), len(name,kind=c_int), memory_resource%c_memory_resource)
end subroutine

subroutine pluto_unregister_resource(name)
    character(len=*), target, intent(in) :: name
    call c_pluto_unregister_resource(c_loc(name), len(name,kind=c_int))
end subroutine

function pluto_host_get_default_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_host_get_default_resource()
end function

function pluto_device_get_default_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_device_get_default_resource()
end function

subroutine pluto_memory_resource_allocate(this, memory, bytes, alignment)
    class(pluto_memory_resource) :: this
    type(c_ptr), intent(out) :: memory
    integer(c_size_t), intent(in) :: bytes
    integer(c_size_t), intent(in), optional :: alignment
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
    if (present(alignment)) then
        call c_pluto_memory_resource_deallocate(this%c_memory_resource, memory, bytes, alignment)
    else
        call c_pluto_memory_resource_deallocate(this%c_memory_resource, memory, bytes, int(0,c_size_t))
    endif
    memory = c_null_ptr
end subroutine

subroutine pluto_memory_pool_resource_release(this)
    class(pluto_memory_resource), intent(in) :: this
    call c_pluto_memory_pool_resource_release(this%c_memory_resource)
end subroutine

subroutine pluto_memory_pool_resource_reserve(this, bytes)
    class(pluto_memory_resource) :: this
    integer(c_size_t), intent(in) :: bytes
    call c_pluto_memory_pool_resource_reserve(this%c_memory_resource, bytes)
end subroutine

function pluto_memory_pool_resource_size(this)
    integer(c_size_t) :: pluto_memory_pool_resource_size
    class(pluto_memory_resource), intent(in) :: this
    pluto_memory_pool_resource_size = c_pluto_memory_pool_resource_size(this%c_memory_resource)
end function

function pluto_memory_pool_resource_capacity(this)
    integer(c_size_t) :: pluto_memory_pool_resource_capacity
    class(pluto_memory_resource), intent(in) :: this
    pluto_memory_pool_resource_capacity = c_pluto_memory_pool_resource_capacity(this%c_memory_resource)
end function

function pluto_host_make_allocator() result(allocator)
    type(pluto_allocator) :: allocator
    allocator%memory_resource = pluto_host_get_default_resource()
end function

function pluto_device_make_allocator() result(allocator)
    type(pluto_allocator) :: allocator
    allocator%memory_resource = pluto_device_get_default_resource()
end function

function pluto_make_allocator_type(this, memory_resource) result(allocator)
    class(pluto_t) :: this
    type(pluto_allocator) :: allocator
    type(pluto_memory_resource) :: memory_resource
    allocator%memory_resource%c_memory_resource = memory_resource%c_memory_resource
end function

function pluto_make_allocator_name(this, memory_resource) result(allocator)
    class(pluto_t) :: this
    type(pluto_allocator) :: allocator
    character(len=*), target, intent(in) :: memory_resource
    allocator%memory_resource%c_memory_resource = &
        & c_pluto_get_registered_resource(c_loc(memory_resource), len(memory_resource,kind=c_int))
end function


subroutine pluto_allocator_allocate_int32_r1(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int32_r2(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int32_r3(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int32_r4(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int64_r1(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int64_r2(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int64_r3(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_int64_r4(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real32_r1(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real32_r2(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real32_r3(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real32_r4(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 4
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real64_r1(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real64_r2(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real64_r3(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_allocate_real64_r4(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = product(shape) * 8
    if (bytes > 0) then
        call this%memory_resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
    else
        array => null()
    endif
end subroutine

subroutine pluto_allocator_deallocate_int32_r1(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_int32_r2(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_int32_r3(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_int32_r4(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1,1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine


subroutine pluto_allocator_deallocate_int64_r1(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_int64_r2(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_int64_r3(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_int64_r4(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1,1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real32_r1(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real32_r2(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real32_r3(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real32_r4(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(1,1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real64_r1(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real64_r2(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real64_r3(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

subroutine pluto_allocator_deallocate_real64_r4(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(1,1,1,1))
        call this%memory_resource%deallocate(mem, bytes)
    endif
    array => null()
end subroutine

function pluto_new_delete_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_new_delete_resource()
end function

function pluto_null_memory_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_null_memory_resource()
end function

function pluto_host_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_host_resource()
end function

function pluto_pinned_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_pinned_resource()
end function

function pluto_device_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_device_resource()
end function

function pluto_managed_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_managed_resource()
end function

function pluto_pinned_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_pinned_pool_resource()
end function

function pluto_host_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_host_pool_resource()
end function

function pluto_device_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_device_pool_resource()
end function

function pluto_managed_pool_resource() result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_managed_pool_resource()
end function

end module
