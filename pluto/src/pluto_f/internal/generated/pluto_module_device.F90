! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_device

use, intrinsic :: iso_c_binding, only : c_int, c_double, c_float, c_int32_t, c_int64_t, c_loc
use pluto_module_allocate_deallocate, only : pluto_allocate, pluto_deallocate
use pluto_module_allocator, only : pluto_allocator
use pluto_module_memory_resource, only : pluto_memory_resource

implicit none
private


! A workaround for NVHPC compiler bug.
! If the compiler is NVHPC we need the THIS argument to be passed explicitly.
! There are otherwise issues with ambiguous generic interfaces between pluto_host_t and pluto_device_t,
! which are resolved by passing the THIS argument explicitly.

#if !defined(__NVCOMPILER)
#define NOPASS , nopass
#define THIS
#define THIS_COMMA
#define CLASS_THIS
#else
#define NOPASS
#define THIS this
#define THIS_COMMA this,
#define CLASS_THIS class(pluto_device_t), intent(in) :: this
#endif

public :: pluto_device_t

type pluto_device_t
   integer :: dummy
contains
    procedure NOPASS :: get_default_resource  => pluto_device_get_default_resource
    procedure NOPASS :: make_allocator        => pluto_device_make_allocator
    procedure NOPASS, private :: pluto_device_set_default_resource_name
    procedure NOPASS, private :: pluto_device_set_default_resource_type
    generic, public :: set_default_resource => pluto_device_set_default_resource_type, &
                                             & pluto_device_set_default_resource_name


    procedure NOPASS, private :: pluto_device_allocate_int32_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_int32_r1_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r1_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int32_r1_bounds, &
        & pluto_device_allocate_int32_r1_shape, &
        & pluto_device_allocate_label_int32_r1_bounds, &
        & pluto_device_allocate_label_int32_r1_shape

    procedure NOPASS, private :: pluto_device_deallocate_int32_r1
    procedure NOPASS, private :: pluto_device_deallocate_label_int32_r1
    generic, public :: deallocate => &
        & pluto_device_deallocate_int32_r1, &
        & pluto_device_deallocate_label_int32_r1
    procedure NOPASS, private :: pluto_device_allocate_int64_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_int64_r1_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r1_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int64_r1_bounds, &
        & pluto_device_allocate_int64_r1_shape, &
        & pluto_device_allocate_label_int64_r1_bounds, &
        & pluto_device_allocate_label_int64_r1_shape

    procedure NOPASS, private :: pluto_device_deallocate_int64_r1
    procedure NOPASS, private :: pluto_device_deallocate_label_int64_r1
    generic, public :: deallocate => &
        & pluto_device_deallocate_int64_r1, &
        & pluto_device_deallocate_label_int64_r1
    procedure NOPASS, private :: pluto_device_allocate_real32_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_real32_r1_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r1_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real32_r1_bounds, &
        & pluto_device_allocate_real32_r1_shape, &
        & pluto_device_allocate_label_real32_r1_bounds, &
        & pluto_device_allocate_label_real32_r1_shape

    procedure NOPASS, private :: pluto_device_deallocate_real32_r1
    procedure NOPASS, private :: pluto_device_deallocate_label_real32_r1
    generic, public :: deallocate => &
        & pluto_device_deallocate_real32_r1, &
        & pluto_device_deallocate_label_real32_r1
    procedure NOPASS, private :: pluto_device_allocate_real64_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_real64_r1_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r1_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r1_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real64_r1_bounds, &
        & pluto_device_allocate_real64_r1_shape, &
        & pluto_device_allocate_label_real64_r1_bounds, &
        & pluto_device_allocate_label_real64_r1_shape

    procedure NOPASS, private :: pluto_device_deallocate_real64_r1
    procedure NOPASS, private :: pluto_device_deallocate_label_real64_r1
    generic, public :: deallocate => &
        & pluto_device_deallocate_real64_r1, &
        & pluto_device_deallocate_label_real64_r1
    procedure NOPASS, private :: pluto_device_allocate_int32_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_int32_r2_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r2_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int32_r2_bounds, &
        & pluto_device_allocate_int32_r2_shape, &
        & pluto_device_allocate_label_int32_r2_bounds, &
        & pluto_device_allocate_label_int32_r2_shape

    procedure NOPASS, private :: pluto_device_deallocate_int32_r2
    procedure NOPASS, private :: pluto_device_deallocate_label_int32_r2
    generic, public :: deallocate => &
        & pluto_device_deallocate_int32_r2, &
        & pluto_device_deallocate_label_int32_r2
    procedure NOPASS, private :: pluto_device_allocate_int64_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_int64_r2_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r2_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int64_r2_bounds, &
        & pluto_device_allocate_int64_r2_shape, &
        & pluto_device_allocate_label_int64_r2_bounds, &
        & pluto_device_allocate_label_int64_r2_shape

    procedure NOPASS, private :: pluto_device_deallocate_int64_r2
    procedure NOPASS, private :: pluto_device_deallocate_label_int64_r2
    generic, public :: deallocate => &
        & pluto_device_deallocate_int64_r2, &
        & pluto_device_deallocate_label_int64_r2
    procedure NOPASS, private :: pluto_device_allocate_real32_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_real32_r2_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r2_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real32_r2_bounds, &
        & pluto_device_allocate_real32_r2_shape, &
        & pluto_device_allocate_label_real32_r2_bounds, &
        & pluto_device_allocate_label_real32_r2_shape

    procedure NOPASS, private :: pluto_device_deallocate_real32_r2
    procedure NOPASS, private :: pluto_device_deallocate_label_real32_r2
    generic, public :: deallocate => &
        & pluto_device_deallocate_real32_r2, &
        & pluto_device_deallocate_label_real32_r2
    procedure NOPASS, private :: pluto_device_allocate_real64_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_real64_r2_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r2_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r2_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real64_r2_bounds, &
        & pluto_device_allocate_real64_r2_shape, &
        & pluto_device_allocate_label_real64_r2_bounds, &
        & pluto_device_allocate_label_real64_r2_shape

    procedure NOPASS, private :: pluto_device_deallocate_real64_r2
    procedure NOPASS, private :: pluto_device_deallocate_label_real64_r2
    generic, public :: deallocate => &
        & pluto_device_deallocate_real64_r2, &
        & pluto_device_deallocate_label_real64_r2
    procedure NOPASS, private :: pluto_device_allocate_int32_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_int32_r3_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r3_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int32_r3_bounds, &
        & pluto_device_allocate_int32_r3_shape, &
        & pluto_device_allocate_label_int32_r3_bounds, &
        & pluto_device_allocate_label_int32_r3_shape

    procedure NOPASS, private :: pluto_device_deallocate_int32_r3
    procedure NOPASS, private :: pluto_device_deallocate_label_int32_r3
    generic, public :: deallocate => &
        & pluto_device_deallocate_int32_r3, &
        & pluto_device_deallocate_label_int32_r3
    procedure NOPASS, private :: pluto_device_allocate_int64_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_int64_r3_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r3_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int64_r3_bounds, &
        & pluto_device_allocate_int64_r3_shape, &
        & pluto_device_allocate_label_int64_r3_bounds, &
        & pluto_device_allocate_label_int64_r3_shape

    procedure NOPASS, private :: pluto_device_deallocate_int64_r3
    procedure NOPASS, private :: pluto_device_deallocate_label_int64_r3
    generic, public :: deallocate => &
        & pluto_device_deallocate_int64_r3, &
        & pluto_device_deallocate_label_int64_r3
    procedure NOPASS, private :: pluto_device_allocate_real32_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_real32_r3_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r3_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real32_r3_bounds, &
        & pluto_device_allocate_real32_r3_shape, &
        & pluto_device_allocate_label_real32_r3_bounds, &
        & pluto_device_allocate_label_real32_r3_shape

    procedure NOPASS, private :: pluto_device_deallocate_real32_r3
    procedure NOPASS, private :: pluto_device_deallocate_label_real32_r3
    generic, public :: deallocate => &
        & pluto_device_deallocate_real32_r3, &
        & pluto_device_deallocate_label_real32_r3
    procedure NOPASS, private :: pluto_device_allocate_real64_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_real64_r3_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r3_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r3_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real64_r3_bounds, &
        & pluto_device_allocate_real64_r3_shape, &
        & pluto_device_allocate_label_real64_r3_bounds, &
        & pluto_device_allocate_label_real64_r3_shape

    procedure NOPASS, private :: pluto_device_deallocate_real64_r3
    procedure NOPASS, private :: pluto_device_deallocate_label_real64_r3
    generic, public :: deallocate => &
        & pluto_device_deallocate_real64_r3, &
        & pluto_device_deallocate_label_real64_r3
    procedure NOPASS, private :: pluto_device_allocate_int32_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_int32_r4_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r4_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int32_r4_bounds, &
        & pluto_device_allocate_int32_r4_shape, &
        & pluto_device_allocate_label_int32_r4_bounds, &
        & pluto_device_allocate_label_int32_r4_shape

    procedure NOPASS, private :: pluto_device_deallocate_int32_r4
    procedure NOPASS, private :: pluto_device_deallocate_label_int32_r4
    generic, public :: deallocate => &
        & pluto_device_deallocate_int32_r4, &
        & pluto_device_deallocate_label_int32_r4
    procedure NOPASS, private :: pluto_device_allocate_int64_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_int64_r4_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r4_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int64_r4_bounds, &
        & pluto_device_allocate_int64_r4_shape, &
        & pluto_device_allocate_label_int64_r4_bounds, &
        & pluto_device_allocate_label_int64_r4_shape

    procedure NOPASS, private :: pluto_device_deallocate_int64_r4
    procedure NOPASS, private :: pluto_device_deallocate_label_int64_r4
    generic, public :: deallocate => &
        & pluto_device_deallocate_int64_r4, &
        & pluto_device_deallocate_label_int64_r4
    procedure NOPASS, private :: pluto_device_allocate_real32_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_real32_r4_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r4_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real32_r4_bounds, &
        & pluto_device_allocate_real32_r4_shape, &
        & pluto_device_allocate_label_real32_r4_bounds, &
        & pluto_device_allocate_label_real32_r4_shape

    procedure NOPASS, private :: pluto_device_deallocate_real32_r4
    procedure NOPASS, private :: pluto_device_deallocate_label_real32_r4
    generic, public :: deallocate => &
        & pluto_device_deallocate_real32_r4, &
        & pluto_device_deallocate_label_real32_r4
    procedure NOPASS, private :: pluto_device_allocate_real64_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_real64_r4_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r4_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r4_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real64_r4_bounds, &
        & pluto_device_allocate_real64_r4_shape, &
        & pluto_device_allocate_label_real64_r4_bounds, &
        & pluto_device_allocate_label_real64_r4_shape

    procedure NOPASS, private :: pluto_device_deallocate_real64_r4
    procedure NOPASS, private :: pluto_device_deallocate_label_real64_r4
    generic, public :: deallocate => &
        & pluto_device_deallocate_real64_r4, &
        & pluto_device_deallocate_label_real64_r4
    procedure NOPASS, private :: pluto_device_allocate_int32_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_int32_r5_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int32_r5_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int32_r5_bounds, &
        & pluto_device_allocate_int32_r5_shape, &
        & pluto_device_allocate_label_int32_r5_bounds, &
        & pluto_device_allocate_label_int32_r5_shape

    procedure NOPASS, private :: pluto_device_deallocate_int32_r5
    procedure NOPASS, private :: pluto_device_deallocate_label_int32_r5
    generic, public :: deallocate => &
        & pluto_device_deallocate_int32_r5, &
        & pluto_device_deallocate_label_int32_r5
    procedure NOPASS, private :: pluto_device_allocate_int64_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_int64_r5_shape
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_int64_r5_shape
    generic, public :: allocate => &
        & pluto_device_allocate_int64_r5_bounds, &
        & pluto_device_allocate_int64_r5_shape, &
        & pluto_device_allocate_label_int64_r5_bounds, &
        & pluto_device_allocate_label_int64_r5_shape

    procedure NOPASS, private :: pluto_device_deallocate_int64_r5
    procedure NOPASS, private :: pluto_device_deallocate_label_int64_r5
    generic, public :: deallocate => &
        & pluto_device_deallocate_int64_r5, &
        & pluto_device_deallocate_label_int64_r5
    procedure NOPASS, private :: pluto_device_allocate_real32_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_real32_r5_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real32_r5_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real32_r5_bounds, &
        & pluto_device_allocate_real32_r5_shape, &
        & pluto_device_allocate_label_real32_r5_bounds, &
        & pluto_device_allocate_label_real32_r5_shape

    procedure NOPASS, private :: pluto_device_deallocate_real32_r5
    procedure NOPASS, private :: pluto_device_deallocate_label_real32_r5
    generic, public :: deallocate => &
        & pluto_device_deallocate_real32_r5, &
        & pluto_device_deallocate_label_real32_r5
    procedure NOPASS, private :: pluto_device_allocate_real64_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_real64_r5_shape
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r5_bounds
    procedure NOPASS, private :: pluto_device_allocate_label_real64_r5_shape
    generic, public :: allocate => &
        & pluto_device_allocate_real64_r5_bounds, &
        & pluto_device_allocate_real64_r5_shape, &
        & pluto_device_allocate_label_real64_r5_bounds, &
        & pluto_device_allocate_label_real64_r5_shape

    procedure NOPASS, private :: pluto_device_deallocate_real64_r5
    procedure NOPASS, private :: pluto_device_deallocate_label_real64_r5
    generic, public :: deallocate => &
        & pluto_device_deallocate_real64_r5, &
        & pluto_device_deallocate_label_real64_r5
end type

contains

function pluto_device_get_default_resource(THIS) result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    CLASS_THIS
    interface
        function c_pluto_device_get_default_resource() result(mr) bind(c)
            use iso_c_binding, only: c_ptr
            type(c_ptr) :: mr 
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_device_get_default_resource()
end function
subroutine pluto_device_set_default_resource_name(THIS_COMMA name)
    CLASS_THIS
    character(len=*), target, intent(in) :: name
    interface
        subroutine c_pluto_device_set_default_resource_name(name, name_size) bind(c)
            use iso_c_binding, only: c_ptr, c_int
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
        end subroutine
    end interface
    call c_pluto_device_set_default_resource_name(c_loc(name), len_trim(name,kind=c_int))
end subroutine

subroutine pluto_device_set_default_resource_type(THIS_COMMA memory_resource)
    CLASS_THIS
    type(pluto_memory_resource), intent(in) :: memory_resource
    interface
        subroutine c_pluto_device_set_default_resource_ptr(mr) bind(c)
            use iso_c_binding, only: c_ptr
            type(c_ptr), value, intent(in) :: mr
        end subroutine
    end interface
    call c_pluto_device_set_default_resource_ptr(memory_resource%c_memory_resource)
end subroutine


function pluto_device_make_allocator(THIS) result(allocator)
    type(pluto_allocator) :: allocator
    CLASS_THIS
    allocator%memory_resource = pluto_device_get_default_resource(THIS)
end function

    subroutine pluto_device_allocate_int32_r1_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r1_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r1_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r1_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int32_r1(THIS_COMMA array)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int32_r1(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r1_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r1_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r1_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r1_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int64_r1(THIS_COMMA array)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int64_r1(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r1_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r1_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r1_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r1_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real32_r1(THIS_COMMA array)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real32_r1(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r1_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r1_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r1_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: lbounds(1)
        integer(c_int), intent(in) :: ubounds(1)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r1_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(1)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real64_r1(THIS_COMMA array)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real64_r1(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r2_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r2_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r2_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r2_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int32_r2(THIS_COMMA array)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int32_r2(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r2_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r2_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r2_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r2_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int64_r2(THIS_COMMA array)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int64_r2(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r2_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r2_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r2_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r2_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real32_r2(THIS_COMMA array)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real32_r2(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r2_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r2_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r2_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: lbounds(2)
        integer(c_int), intent(in) :: ubounds(2)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r2_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(2)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real64_r2(THIS_COMMA array)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real64_r2(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r3_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r3_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r3_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r3_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int32_r3(THIS_COMMA array)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int32_r3(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r3_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r3_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r3_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r3_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int64_r3(THIS_COMMA array)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int64_r3(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r3_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r3_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r3_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r3_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real32_r3(THIS_COMMA array)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real32_r3(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r3_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r3_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r3_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: lbounds(3)
        integer(c_int), intent(in) :: ubounds(3)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r3_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(3)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real64_r3(THIS_COMMA array)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real64_r3(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r4_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r4_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r4_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r4_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int32_r4(THIS_COMMA array)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int32_r4(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r4_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r4_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r4_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r4_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int64_r4(THIS_COMMA array)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int64_r4(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r4_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r4_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r4_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r4_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real32_r4(THIS_COMMA array)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real32_r4(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r4_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r4_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r4_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: lbounds(4)
        integer(c_int), intent(in) :: ubounds(4)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r4_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(4)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real64_r4(THIS_COMMA array)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real64_r4(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r5_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int32_r5_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r5_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int32_r5_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int32_r5(THIS_COMMA array)
        CLASS_THIS
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int32_r5(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r5_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_int64_r5_shape(THIS_COMMA array, shape)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r5_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_int64_r5_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_int64_r5(THIS_COMMA array)
        CLASS_THIS
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_int64_r5(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r5_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real32_r5_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r5_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real32_r5_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real32_r5(THIS_COMMA array)
        CLASS_THIS
        real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real32_r5(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r5_bounds(THIS_COMMA array, lbounds, ubounds)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_real64_r5_shape(THIS_COMMA array, shape)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r5_bounds(THIS_COMMA label, array, lbounds, ubounds)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: lbounds(5)
        integer(c_int), intent(in) :: ubounds(5)
        call pluto_allocate(label, array, lbounds, ubounds, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_allocate_label_real64_r5_shape(THIS_COMMA label, array, shape)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(c_int), intent(in) :: shape(5)
        call pluto_allocate(label, array, shape, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_real64_r5(THIS_COMMA array)
        CLASS_THIS
        real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(array, pluto_device_get_default_resource(THIS))
    end subroutine
    subroutine pluto_device_deallocate_label_real64_r5(THIS_COMMA label, array)
        CLASS_THIS
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
        call pluto_deallocate(label, array, pluto_device_get_default_resource(THIS))
    end subroutine

end module
