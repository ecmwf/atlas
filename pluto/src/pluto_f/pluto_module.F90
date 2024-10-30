module pluto_module

use iso_c_binding, only: c_loc, c_ptr, c_int, c_size_t, c_null_ptr, c_double, c_float, c_int32_t, c_int64_t, c_f_pointer
implicit none
private

public :: pluto, pluto_memory_resource, pluto_allocator

interface
    subroutine c_pluto_host_set_default_resource(name, name_size) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end subroutine
    subroutine c_pluto_device_set_default_resource(name, name_size) bind(c)
        use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: name
        integer(c_int), value, intent(in) :: name_size
    end subroutine
    function c_pluto_host_get_default_resource() result(memory_resource) bind(c)
        use iso_c_binding, only: c_ptr
        type(c_ptr) :: memory_resource
    end function
    function c_pluto_device_get_default_resource() result(memory_resource) bind(c)
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
    subroutine c_pluto_scope_push() bind(c)
    end subroutine
    subroutine c_pluto_scope_pop() bind(c)
    end subroutine
end interface

type pluto_memory_resource
    type(c_ptr) :: c_memory_resource
contains
    procedure :: allocate => pluto_memory_resource_allocate
    procedure :: deallocate => pluto_memory_resource_deallocate
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
    procedure :: pluto_allocator_allocate_mold_int32_r1
    procedure :: pluto_allocator_allocate_mold_int32_r2
    procedure :: pluto_allocator_allocate_mold_int32_r3
    procedure :: pluto_allocator_allocate_mold_int32_r4
    procedure :: pluto_allocator_allocate_mold_int64_r1
    procedure :: pluto_allocator_allocate_mold_int64_r2
    procedure :: pluto_allocator_allocate_mold_int64_r3
    procedure :: pluto_allocator_allocate_mold_int64_r4
    procedure :: pluto_allocator_allocate_mold_real32_r1
    procedure :: pluto_allocator_allocate_mold_real32_r2
    procedure :: pluto_allocator_allocate_mold_real32_r3
    procedure :: pluto_allocator_allocate_mold_real32_r4
    procedure :: pluto_allocator_allocate_mold_real64_r1
    procedure :: pluto_allocator_allocate_mold_real64_r2
    procedure :: pluto_allocator_allocate_mold_real64_r3
    procedure :: pluto_allocator_allocate_mold_real64_r4

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
        & pluto_allocator_allocate_real64_r4, &
        & pluto_allocator_allocate_mold_int32_r1, &
        & pluto_allocator_allocate_mold_int32_r2, &
        & pluto_allocator_allocate_mold_int32_r3, &
        & pluto_allocator_allocate_mold_int32_r4, &
        & pluto_allocator_allocate_mold_int64_r1, &
        & pluto_allocator_allocate_mold_int64_r2, &
        & pluto_allocator_allocate_mold_int64_r3, &
        & pluto_allocator_allocate_mold_int64_r4, &
        & pluto_allocator_allocate_mold_real32_r1, &
        & pluto_allocator_allocate_mold_real32_r2, &
        & pluto_allocator_allocate_mold_real32_r3, &
        & pluto_allocator_allocate_mold_real32_r4, &
        & pluto_allocator_allocate_mold_real64_r1, &
        & pluto_allocator_allocate_mold_real64_r2, &
        & pluto_allocator_allocate_mold_real64_r3, &
        & pluto_allocator_allocate_mold_real64_r4


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
contains
    procedure, nopass :: set_default_resource  => pluto_host_set_default_resource
    procedure, nopass :: get_default_resource  => pluto_host_get_default_resource
    procedure, nopass :: get_default_allocator => pluto_host_get_default_allocator
end type

type pluto_device_t
contains
    procedure, nopass :: set_default_resource  => pluto_device_set_default_resource
    procedure, nopass :: get_default_resource  => pluto_device_get_default_resource
    procedure, nopass :: get_default_allocator => pluto_device_get_default_allocator
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
end type

type(pluto_t) :: pluto

contains

subroutine pluto_host_set_default_resource(name)
    character(len=*), target, intent(in) :: name
    call c_pluto_host_set_default_resource(c_loc(name), len(name,kind=c_int))
end subroutine

subroutine pluto_device_set_default_resource(name)
    character(len=*), target, intent(in) :: name
    call c_pluto_device_set_default_resource(c_loc(name), len(name,kind=c_int))
end subroutine

subroutine pluto_scope_push()
    call c_pluto_scope_push()
end subroutine

subroutine pluto_scope_pop()
    call c_pluto_scope_pop()
end subroutine

subroutine pluto_host_get_default_resource(memory_resource)
    type(pluto_memory_resource), intent(out) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_host_get_default_resource()
end subroutine

subroutine pluto_device_get_default_resource(memory_resource)
    type(pluto_memory_resource), intent(out) :: memory_resource
    memory_resource%c_memory_resource = c_pluto_device_get_default_resource()
end subroutine

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

subroutine pluto_host_get_default_allocator(allocator)
    type(pluto_allocator) :: allocator
    call pluto_host_get_default_resource(allocator%memory_resource)
end subroutine

subroutine pluto_device_get_default_allocator(allocator)
    type(pluto_allocator) :: allocator
    call pluto_device_get_default_resource(allocator%memory_resource)
end subroutine

subroutine pluto_allocator_allocate_int32_r1(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(out) :: array(:)
    integer(c_int), intent(in) :: shape(1)
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
    integer(c_int32_t), pointer, intent(out) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
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
    integer(c_int32_t), pointer, intent(out) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
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
    integer(c_int32_t), pointer, intent(out) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
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
    integer(c_int64_t), pointer, intent(out) :: array(:)
    integer(c_int), intent(in) :: shape(1)
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
    integer(c_int64_t), pointer, intent(out) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
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
    integer(c_int64_t), pointer, intent(out) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
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
    integer(c_int64_t), pointer, intent(out) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
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
    real(c_float), pointer, intent(out) :: array(:)
    integer(c_int), intent(in) :: shape(1)
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
    real(c_float), pointer, intent(out) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
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
    real(c_float), pointer, intent(out) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
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
    real(c_float), pointer, intent(out) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
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
    real(c_double), pointer, intent(out) :: array(:)
    integer(c_int), intent(in) :: shape(1)
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
    real(c_double), pointer, intent(out) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
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
    real(c_double), pointer, intent(out) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
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
    real(c_double), pointer, intent(out) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
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


subroutine pluto_allocator_allocate_mold_int32_r1(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(out) :: array(:)
    integer(c_int32_t), intent(in) :: mold(:)
    call pluto_allocator_allocate_int32_r1(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int32_r2(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(out) :: array(:,:)
    integer(c_int32_t), intent(in) :: mold(:,:)
    call pluto_allocator_allocate_int32_r2(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int32_r3(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(out) :: array(:,:,:)
    integer(c_int32_t), intent(in) :: mold(:,:,:)
    call pluto_allocator_allocate_int32_r3(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int32_r4(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(out) :: array(:,:,:,:)
    integer(c_int32_t), intent(in) :: mold(:,:,:,:)
    call pluto_allocator_allocate_int32_r4(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int64_r1(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(out) :: array(:)
    integer(c_int64_t), intent(in) :: mold(:)
    call pluto_allocator_allocate_int64_r1(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int64_r2(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(out) :: array(:,:)
    integer(c_int64_t), intent(in) :: mold(:,:)
    call pluto_allocator_allocate_int64_r2(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int64_r3(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(out) :: array(:,:,:)
    integer(c_int64_t), intent(in) :: mold(:,:,:)
    call pluto_allocator_allocate_int64_r3(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_int64_r4(this, array, mold)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(out) :: array(:,:,:,:)
    integer(c_int64_t), intent(in) :: mold(:,:,:,:)
    call pluto_allocator_allocate_int64_r4(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real32_r1(this, array, mold)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(out) :: array(:)
    real(c_float), intent(in) :: mold(:)
    call pluto_allocator_allocate_real32_r1(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real32_r2(this, array, mold)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(out) :: array(:,:)
    real(c_float), intent(in) :: mold(:,:)
    call pluto_allocator_allocate_real32_r2(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real32_r3(this, array, mold)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(out) :: array(:,:,:)
    real(c_float), intent(in) :: mold(:,:,:)
    call pluto_allocator_allocate_real32_r3(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real32_r4(this, array, mold)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(out) :: array(:,:,:,:)
    real(c_float), intent(in) :: mold(:,:,:,:)
    call pluto_allocator_allocate_real32_r4(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real64_r1(this, array, mold)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(out) :: array(:)
    real(c_double), intent(in) :: mold(:)
    call pluto_allocator_allocate_real64_r1(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real64_r2(this, array, mold)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(out) :: array(:,:)
    real(c_double), intent(in) :: mold(:,:)
    call pluto_allocator_allocate_real64_r2(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real64_r3(this, array, mold)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(out) :: array(:,:,:)
    real(c_double), intent(in) :: mold(:,:,:)
    call pluto_allocator_allocate_real64_r3(this, array, shape(mold))
end subroutine

subroutine pluto_allocator_allocate_mold_real64_r4(this, array, mold)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(out) :: array(:,:,:,:)
    real(c_double), intent(in) :: mold(:,:,:,:)
    call pluto_allocator_allocate_real64_r4(this, array, shape(mold))
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

end module
