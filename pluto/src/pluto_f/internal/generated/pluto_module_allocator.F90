! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_allocator

use, intrinsic :: iso_c_binding, only : c_loc, c_ptr, c_int, c_size_t, c_null_ptr, c_double, c_float, &
                                      & c_int32_t, c_int64_t, c_f_pointer, c_associated
use pluto_module_abort, only : &
    pluto_abort
use pluto_module_memory_resource, only : &
    pluto_memory_resource, &
    pluto_get_registered_resource
use pluto_module_allocate_deallocate, only : &
    pluto_allocate, &
    pluto_deallocate

implicit none
private


public :: pluto_allocator
public :: pluto_make_allocator

type pluto_allocator
    type(pluto_memory_resource) :: memory_resource
contains

    procedure, private :: pluto_allocator_allocate_int32_r1_shape
    procedure, private :: pluto_allocator_allocate_int32_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r1_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r1_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int32_r1_shape, &
        & pluto_allocator_allocate_int32_r1_bounds, &
        & pluto_allocator_allocate_label_int32_r1_shape, &
        & pluto_allocator_allocate_label_int32_r1_bounds

    procedure, private :: pluto_allocator_deallocate_int32_r1
    procedure, private :: pluto_allocator_deallocate_label_int32_r1

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int32_r1, &
        & pluto_allocator_deallocate_label_int32_r1
    procedure, private :: pluto_allocator_allocate_int64_r1_shape
    procedure, private :: pluto_allocator_allocate_int64_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r1_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r1_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int64_r1_shape, &
        & pluto_allocator_allocate_int64_r1_bounds, &
        & pluto_allocator_allocate_label_int64_r1_shape, &
        & pluto_allocator_allocate_label_int64_r1_bounds

    procedure, private :: pluto_allocator_deallocate_int64_r1
    procedure, private :: pluto_allocator_deallocate_label_int64_r1

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int64_r1, &
        & pluto_allocator_deallocate_label_int64_r1
    procedure, private :: pluto_allocator_allocate_real32_r1_shape
    procedure, private :: pluto_allocator_allocate_real32_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r1_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r1_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real32_r1_shape, &
        & pluto_allocator_allocate_real32_r1_bounds, &
        & pluto_allocator_allocate_label_real32_r1_shape, &
        & pluto_allocator_allocate_label_real32_r1_bounds

    procedure, private :: pluto_allocator_deallocate_real32_r1
    procedure, private :: pluto_allocator_deallocate_label_real32_r1

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real32_r1, &
        & pluto_allocator_deallocate_label_real32_r1
    procedure, private :: pluto_allocator_allocate_real64_r1_shape
    procedure, private :: pluto_allocator_allocate_real64_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r1_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r1_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real64_r1_shape, &
        & pluto_allocator_allocate_real64_r1_bounds, &
        & pluto_allocator_allocate_label_real64_r1_shape, &
        & pluto_allocator_allocate_label_real64_r1_bounds

    procedure, private :: pluto_allocator_deallocate_real64_r1
    procedure, private :: pluto_allocator_deallocate_label_real64_r1

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real64_r1, &
        & pluto_allocator_deallocate_label_real64_r1
    procedure, private :: pluto_allocator_allocate_int32_r2_shape
    procedure, private :: pluto_allocator_allocate_int32_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r2_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r2_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int32_r2_shape, &
        & pluto_allocator_allocate_int32_r2_bounds, &
        & pluto_allocator_allocate_label_int32_r2_shape, &
        & pluto_allocator_allocate_label_int32_r2_bounds

    procedure, private :: pluto_allocator_deallocate_int32_r2
    procedure, private :: pluto_allocator_deallocate_label_int32_r2

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int32_r2, &
        & pluto_allocator_deallocate_label_int32_r2
    procedure, private :: pluto_allocator_allocate_int64_r2_shape
    procedure, private :: pluto_allocator_allocate_int64_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r2_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r2_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int64_r2_shape, &
        & pluto_allocator_allocate_int64_r2_bounds, &
        & pluto_allocator_allocate_label_int64_r2_shape, &
        & pluto_allocator_allocate_label_int64_r2_bounds

    procedure, private :: pluto_allocator_deallocate_int64_r2
    procedure, private :: pluto_allocator_deallocate_label_int64_r2

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int64_r2, &
        & pluto_allocator_deallocate_label_int64_r2
    procedure, private :: pluto_allocator_allocate_real32_r2_shape
    procedure, private :: pluto_allocator_allocate_real32_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r2_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r2_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real32_r2_shape, &
        & pluto_allocator_allocate_real32_r2_bounds, &
        & pluto_allocator_allocate_label_real32_r2_shape, &
        & pluto_allocator_allocate_label_real32_r2_bounds

    procedure, private :: pluto_allocator_deallocate_real32_r2
    procedure, private :: pluto_allocator_deallocate_label_real32_r2

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real32_r2, &
        & pluto_allocator_deallocate_label_real32_r2
    procedure, private :: pluto_allocator_allocate_real64_r2_shape
    procedure, private :: pluto_allocator_allocate_real64_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r2_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r2_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real64_r2_shape, &
        & pluto_allocator_allocate_real64_r2_bounds, &
        & pluto_allocator_allocate_label_real64_r2_shape, &
        & pluto_allocator_allocate_label_real64_r2_bounds

    procedure, private :: pluto_allocator_deallocate_real64_r2
    procedure, private :: pluto_allocator_deallocate_label_real64_r2

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real64_r2, &
        & pluto_allocator_deallocate_label_real64_r2
    procedure, private :: pluto_allocator_allocate_int32_r3_shape
    procedure, private :: pluto_allocator_allocate_int32_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r3_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r3_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int32_r3_shape, &
        & pluto_allocator_allocate_int32_r3_bounds, &
        & pluto_allocator_allocate_label_int32_r3_shape, &
        & pluto_allocator_allocate_label_int32_r3_bounds

    procedure, private :: pluto_allocator_deallocate_int32_r3
    procedure, private :: pluto_allocator_deallocate_label_int32_r3

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int32_r3, &
        & pluto_allocator_deallocate_label_int32_r3
    procedure, private :: pluto_allocator_allocate_int64_r3_shape
    procedure, private :: pluto_allocator_allocate_int64_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r3_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r3_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int64_r3_shape, &
        & pluto_allocator_allocate_int64_r3_bounds, &
        & pluto_allocator_allocate_label_int64_r3_shape, &
        & pluto_allocator_allocate_label_int64_r3_bounds

    procedure, private :: pluto_allocator_deallocate_int64_r3
    procedure, private :: pluto_allocator_deallocate_label_int64_r3

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int64_r3, &
        & pluto_allocator_deallocate_label_int64_r3
    procedure, private :: pluto_allocator_allocate_real32_r3_shape
    procedure, private :: pluto_allocator_allocate_real32_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r3_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r3_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real32_r3_shape, &
        & pluto_allocator_allocate_real32_r3_bounds, &
        & pluto_allocator_allocate_label_real32_r3_shape, &
        & pluto_allocator_allocate_label_real32_r3_bounds

    procedure, private :: pluto_allocator_deallocate_real32_r3
    procedure, private :: pluto_allocator_deallocate_label_real32_r3

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real32_r3, &
        & pluto_allocator_deallocate_label_real32_r3
    procedure, private :: pluto_allocator_allocate_real64_r3_shape
    procedure, private :: pluto_allocator_allocate_real64_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r3_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r3_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real64_r3_shape, &
        & pluto_allocator_allocate_real64_r3_bounds, &
        & pluto_allocator_allocate_label_real64_r3_shape, &
        & pluto_allocator_allocate_label_real64_r3_bounds

    procedure, private :: pluto_allocator_deallocate_real64_r3
    procedure, private :: pluto_allocator_deallocate_label_real64_r3

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real64_r3, &
        & pluto_allocator_deallocate_label_real64_r3
    procedure, private :: pluto_allocator_allocate_int32_r4_shape
    procedure, private :: pluto_allocator_allocate_int32_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r4_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r4_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int32_r4_shape, &
        & pluto_allocator_allocate_int32_r4_bounds, &
        & pluto_allocator_allocate_label_int32_r4_shape, &
        & pluto_allocator_allocate_label_int32_r4_bounds

    procedure, private :: pluto_allocator_deallocate_int32_r4
    procedure, private :: pluto_allocator_deallocate_label_int32_r4

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int32_r4, &
        & pluto_allocator_deallocate_label_int32_r4
    procedure, private :: pluto_allocator_allocate_int64_r4_shape
    procedure, private :: pluto_allocator_allocate_int64_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r4_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r4_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int64_r4_shape, &
        & pluto_allocator_allocate_int64_r4_bounds, &
        & pluto_allocator_allocate_label_int64_r4_shape, &
        & pluto_allocator_allocate_label_int64_r4_bounds

    procedure, private :: pluto_allocator_deallocate_int64_r4
    procedure, private :: pluto_allocator_deallocate_label_int64_r4

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int64_r4, &
        & pluto_allocator_deallocate_label_int64_r4
    procedure, private :: pluto_allocator_allocate_real32_r4_shape
    procedure, private :: pluto_allocator_allocate_real32_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r4_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r4_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real32_r4_shape, &
        & pluto_allocator_allocate_real32_r4_bounds, &
        & pluto_allocator_allocate_label_real32_r4_shape, &
        & pluto_allocator_allocate_label_real32_r4_bounds

    procedure, private :: pluto_allocator_deallocate_real32_r4
    procedure, private :: pluto_allocator_deallocate_label_real32_r4

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real32_r4, &
        & pluto_allocator_deallocate_label_real32_r4
    procedure, private :: pluto_allocator_allocate_real64_r4_shape
    procedure, private :: pluto_allocator_allocate_real64_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r4_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r4_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real64_r4_shape, &
        & pluto_allocator_allocate_real64_r4_bounds, &
        & pluto_allocator_allocate_label_real64_r4_shape, &
        & pluto_allocator_allocate_label_real64_r4_bounds

    procedure, private :: pluto_allocator_deallocate_real64_r4
    procedure, private :: pluto_allocator_deallocate_label_real64_r4

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real64_r4, &
        & pluto_allocator_deallocate_label_real64_r4
    procedure, private :: pluto_allocator_allocate_int32_r5_shape
    procedure, private :: pluto_allocator_allocate_int32_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r5_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r5_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int32_r5_shape, &
        & pluto_allocator_allocate_int32_r5_bounds, &
        & pluto_allocator_allocate_label_int32_r5_shape, &
        & pluto_allocator_allocate_label_int32_r5_bounds

    procedure, private :: pluto_allocator_deallocate_int32_r5
    procedure, private :: pluto_allocator_deallocate_label_int32_r5

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int32_r5, &
        & pluto_allocator_deallocate_label_int32_r5
    procedure, private :: pluto_allocator_allocate_int64_r5_shape
    procedure, private :: pluto_allocator_allocate_int64_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r5_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r5_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_int64_r5_shape, &
        & pluto_allocator_allocate_int64_r5_bounds, &
        & pluto_allocator_allocate_label_int64_r5_shape, &
        & pluto_allocator_allocate_label_int64_r5_bounds

    procedure, private :: pluto_allocator_deallocate_int64_r5
    procedure, private :: pluto_allocator_deallocate_label_int64_r5

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_int64_r5, &
        & pluto_allocator_deallocate_label_int64_r5
    procedure, private :: pluto_allocator_allocate_real32_r5_shape
    procedure, private :: pluto_allocator_allocate_real32_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r5_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r5_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real32_r5_shape, &
        & pluto_allocator_allocate_real32_r5_bounds, &
        & pluto_allocator_allocate_label_real32_r5_shape, &
        & pluto_allocator_allocate_label_real32_r5_bounds

    procedure, private :: pluto_allocator_deallocate_real32_r5
    procedure, private :: pluto_allocator_deallocate_label_real32_r5

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real32_r5, &
        & pluto_allocator_deallocate_label_real32_r5
    procedure, private :: pluto_allocator_allocate_real64_r5_shape
    procedure, private :: pluto_allocator_allocate_real64_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r5_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r5_bounds

    generic, public :: allocate => &
        & pluto_allocator_allocate_real64_r5_shape, &
        & pluto_allocator_allocate_real64_r5_bounds, &
        & pluto_allocator_allocate_label_real64_r5_shape, &
        & pluto_allocator_allocate_label_real64_r5_bounds

    procedure, private :: pluto_allocator_deallocate_real64_r5
    procedure, private :: pluto_allocator_deallocate_label_real64_r5

    generic, public :: deallocate => &
        & pluto_allocator_deallocate_real64_r5, &
        & pluto_allocator_deallocate_label_real64_r5
end type

interface pluto_make_allocator
    module procedure pluto_make_allocator_type
    module procedure pluto_make_allocator_name
end interface

contains

subroutine assert_allocator_is_setup(allocator)
    type(pluto_allocator) :: allocator
    if (.not. c_associated(allocator%memory_resource%c_memory_resource)) then
        call pluto_abort("pluto_allocator has not been assigned or setup properly.&
          & Please ensure that the allocator is created by pluto%make_allocator or&
          & pluto%{host,device}%make_allocator.")
    endif
end subroutine

function pluto_make_allocator_type(resource) result(allocator)
    type(pluto_allocator) :: allocator
    type(pluto_memory_resource) :: resource
    allocator%memory_resource%c_memory_resource = resource%c_memory_resource
end function

function pluto_make_allocator_name(resource) result(allocator)
    type(pluto_allocator) :: allocator
    character(len=*), target, intent(in) :: resource
    allocator%memory_resource = pluto_get_registered_resource(resource)
end function


subroutine pluto_allocator_allocate_int32_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int32_r1(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int32_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int64_r1(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int64_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real32_r1(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real32_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real64_r1(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real64_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int32_r2(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int32_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int64_r2(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int64_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real32_r2(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real32_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real64_r2(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real64_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int32_r3(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int32_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int64_r3(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int64_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real32_r3(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real32_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real64_r3(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real64_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int32_r4(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int32_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int64_r4(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int64_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real32_r4(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real32_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real64_r4(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real64_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int32_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int32_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int32_r5(this, array)
    class(pluto_allocator) :: this
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int32_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_int64_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_int64_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_int64_r5(this, array)
    class(pluto_allocator) :: this
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_int64_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real32_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real32_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real32_r5(this, array)
    class(pluto_allocator) :: this
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real32_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_real64_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end subroutine
subroutine pluto_allocator_allocate_label_real64_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_real64_r5(this, array)
    class(pluto_allocator) :: this
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end subroutine
subroutine pluto_allocator_deallocate_label_real64_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end subroutine

end module
