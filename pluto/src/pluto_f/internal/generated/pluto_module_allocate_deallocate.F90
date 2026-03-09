! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_allocate_deallocate
! This is a separate module which is used by pluto_module to implement the allocation and deallocation procedures
! which use the C memory resource interface.

use, intrinsic :: iso_c_binding, only : c_loc, c_ptr, c_int, c_size_t, c_null_ptr, c_double, c_float, &
                                      & c_int32_t, c_int64_t, c_f_pointer, c_associated, c_f_pointer
use pluto_module_memory_resource, only : pluto_memory_resource, pluto_get_label, pluto_set_label
implicit none
private


public :: pluto_allocate, pluto_deallocate

interface pluto_allocate
    module procedure pluto_allocate_int32_r1_bounds
    module procedure pluto_allocate_int32_r1_shape
    module procedure pluto_allocate_label_int32_r1_bounds
    module procedure pluto_allocate_label_int32_r1_shape
    module procedure pluto_allocate_int64_r1_bounds
    module procedure pluto_allocate_int64_r1_shape
    module procedure pluto_allocate_label_int64_r1_bounds
    module procedure pluto_allocate_label_int64_r1_shape
    module procedure pluto_allocate_real32_r1_bounds
    module procedure pluto_allocate_real32_r1_shape
    module procedure pluto_allocate_label_real32_r1_bounds
    module procedure pluto_allocate_label_real32_r1_shape
    module procedure pluto_allocate_real64_r1_bounds
    module procedure pluto_allocate_real64_r1_shape
    module procedure pluto_allocate_label_real64_r1_bounds
    module procedure pluto_allocate_label_real64_r1_shape
    module procedure pluto_allocate_int32_r2_bounds
    module procedure pluto_allocate_int32_r2_shape
    module procedure pluto_allocate_label_int32_r2_bounds
    module procedure pluto_allocate_label_int32_r2_shape
    module procedure pluto_allocate_int64_r2_bounds
    module procedure pluto_allocate_int64_r2_shape
    module procedure pluto_allocate_label_int64_r2_bounds
    module procedure pluto_allocate_label_int64_r2_shape
    module procedure pluto_allocate_real32_r2_bounds
    module procedure pluto_allocate_real32_r2_shape
    module procedure pluto_allocate_label_real32_r2_bounds
    module procedure pluto_allocate_label_real32_r2_shape
    module procedure pluto_allocate_real64_r2_bounds
    module procedure pluto_allocate_real64_r2_shape
    module procedure pluto_allocate_label_real64_r2_bounds
    module procedure pluto_allocate_label_real64_r2_shape
    module procedure pluto_allocate_int32_r3_bounds
    module procedure pluto_allocate_int32_r3_shape
    module procedure pluto_allocate_label_int32_r3_bounds
    module procedure pluto_allocate_label_int32_r3_shape
    module procedure pluto_allocate_int64_r3_bounds
    module procedure pluto_allocate_int64_r3_shape
    module procedure pluto_allocate_label_int64_r3_bounds
    module procedure pluto_allocate_label_int64_r3_shape
    module procedure pluto_allocate_real32_r3_bounds
    module procedure pluto_allocate_real32_r3_shape
    module procedure pluto_allocate_label_real32_r3_bounds
    module procedure pluto_allocate_label_real32_r3_shape
    module procedure pluto_allocate_real64_r3_bounds
    module procedure pluto_allocate_real64_r3_shape
    module procedure pluto_allocate_label_real64_r3_bounds
    module procedure pluto_allocate_label_real64_r3_shape
    module procedure pluto_allocate_int32_r4_bounds
    module procedure pluto_allocate_int32_r4_shape
    module procedure pluto_allocate_label_int32_r4_bounds
    module procedure pluto_allocate_label_int32_r4_shape
    module procedure pluto_allocate_int64_r4_bounds
    module procedure pluto_allocate_int64_r4_shape
    module procedure pluto_allocate_label_int64_r4_bounds
    module procedure pluto_allocate_label_int64_r4_shape
    module procedure pluto_allocate_real32_r4_bounds
    module procedure pluto_allocate_real32_r4_shape
    module procedure pluto_allocate_label_real32_r4_bounds
    module procedure pluto_allocate_label_real32_r4_shape
    module procedure pluto_allocate_real64_r4_bounds
    module procedure pluto_allocate_real64_r4_shape
    module procedure pluto_allocate_label_real64_r4_bounds
    module procedure pluto_allocate_label_real64_r4_shape
    module procedure pluto_allocate_int32_r5_bounds
    module procedure pluto_allocate_int32_r5_shape
    module procedure pluto_allocate_label_int32_r5_bounds
    module procedure pluto_allocate_label_int32_r5_shape
    module procedure pluto_allocate_int64_r5_bounds
    module procedure pluto_allocate_int64_r5_shape
    module procedure pluto_allocate_label_int64_r5_bounds
    module procedure pluto_allocate_label_int64_r5_shape
    module procedure pluto_allocate_real32_r5_bounds
    module procedure pluto_allocate_real32_r5_shape
    module procedure pluto_allocate_label_real32_r5_bounds
    module procedure pluto_allocate_label_real32_r5_shape
    module procedure pluto_allocate_real64_r5_bounds
    module procedure pluto_allocate_real64_r5_shape
    module procedure pluto_allocate_label_real64_r5_bounds
    module procedure pluto_allocate_label_real64_r5_shape
end interface

interface pluto_deallocate
    module procedure pluto_deallocate_int32_r1
    module procedure pluto_deallocate_label_int32_r1
    module procedure pluto_deallocate_int64_r1
    module procedure pluto_deallocate_label_int64_r1
    module procedure pluto_deallocate_real32_r1
    module procedure pluto_deallocate_label_real32_r1
    module procedure pluto_deallocate_real64_r1
    module procedure pluto_deallocate_label_real64_r1
    module procedure pluto_deallocate_int32_r2
    module procedure pluto_deallocate_label_int32_r2
    module procedure pluto_deallocate_int64_r2
    module procedure pluto_deallocate_label_int64_r2
    module procedure pluto_deallocate_real32_r2
    module procedure pluto_deallocate_label_real32_r2
    module procedure pluto_deallocate_real64_r2
    module procedure pluto_deallocate_label_real64_r2
    module procedure pluto_deallocate_int32_r3
    module procedure pluto_deallocate_label_int32_r3
    module procedure pluto_deallocate_int64_r3
    module procedure pluto_deallocate_label_int64_r3
    module procedure pluto_deallocate_real32_r3
    module procedure pluto_deallocate_label_real32_r3
    module procedure pluto_deallocate_real64_r3
    module procedure pluto_deallocate_label_real64_r3
    module procedure pluto_deallocate_int32_r4
    module procedure pluto_deallocate_label_int32_r4
    module procedure pluto_deallocate_int64_r4
    module procedure pluto_deallocate_label_int64_r4
    module procedure pluto_deallocate_real32_r4
    module procedure pluto_deallocate_label_real32_r4
    module procedure pluto_deallocate_real64_r4
    module procedure pluto_deallocate_label_real64_r4
    module procedure pluto_deallocate_int32_r5
    module procedure pluto_deallocate_label_int32_r5
    module procedure pluto_deallocate_int64_r5
    module procedure pluto_deallocate_label_int64_r5
    module procedure pluto_deallocate_real32_r5
    module procedure pluto_deallocate_label_real32_r5
    module procedure pluto_deallocate_real64_r5
    module procedure pluto_deallocate_label_real64_r5
end interface

contains

function array_size(lbounds, ubounds) result(num)
    implicit none
    integer(c_int), intent(in) :: lbounds(:)
    integer(c_int), intent(in) :: ubounds(:)
    integer(c_int) :: num
    integer(c_int) :: shape(size(lbounds))
    shape = ubounds - lbounds + 1
    num = product(shape)
end function


subroutine pluto_allocate_int32_r1_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(1)
    integer(c_int32_t), pointer :: array_standard_bounds(:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):) => array_standard_bounds(:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int32_r1_shape(array, shape, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_int32_r1_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r1_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int32_r1_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int32_r1(array, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int32_r1(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r1(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int64_r1_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(1)
    integer(c_int64_t), pointer :: array_standard_bounds(:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):) => array_standard_bounds(:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int64_r1_shape(array, shape, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_int64_r1_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r1_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int64_r1_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int64_r1(array, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int64_r1(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r1(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real32_r1_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(1)
    real(c_float), pointer :: array_standard_bounds(:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):) => array_standard_bounds(:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real32_r1_shape(array, shape, resource)
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_real32_r1_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r1_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real32_r1_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real32_r1(array, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real32_r1(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r1(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real64_r1_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(1)
    real(c_double), pointer :: array_standard_bounds(:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):) => array_standard_bounds(:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real64_r1_shape(array, shape, resource)
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_real64_r1_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r1_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real64_r1_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real64_r1(array, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real64_r1(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r1(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int32_r2_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(2)
    integer(c_int32_t), pointer :: array_standard_bounds(:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):) => array_standard_bounds(:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int32_r2_shape(array, shape, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_int32_r2_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r2_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int32_r2_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int32_r2(array, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int32_r2(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r2(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int64_r2_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(2)
    integer(c_int64_t), pointer :: array_standard_bounds(:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):) => array_standard_bounds(:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int64_r2_shape(array, shape, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_int64_r2_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r2_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int64_r2_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int64_r2(array, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int64_r2(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r2(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real32_r2_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(2)
    real(c_float), pointer :: array_standard_bounds(:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):) => array_standard_bounds(:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real32_r2_shape(array, shape, resource)
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_real32_r2_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r2_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real32_r2_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real32_r2(array, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real32_r2(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r2(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real64_r2_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(2)
    real(c_double), pointer :: array_standard_bounds(:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):) => array_standard_bounds(:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real64_r2_shape(array, shape, resource)
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_real64_r2_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r2_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real64_r2_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real64_r2(array, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real64_r2(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r2(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int32_r3_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(3)
    integer(c_int32_t), pointer :: array_standard_bounds(:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array_standard_bounds(:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int32_r3_shape(array, shape, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_int32_r3_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r3_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int32_r3_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int32_r3(array, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int32_r3(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r3(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int64_r3_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(3)
    integer(c_int64_t), pointer :: array_standard_bounds(:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array_standard_bounds(:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int64_r3_shape(array, shape, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_int64_r3_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r3_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int64_r3_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int64_r3(array, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int64_r3(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r3(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real32_r3_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(3)
    real(c_float), pointer :: array_standard_bounds(:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array_standard_bounds(:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real32_r3_shape(array, shape, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_real32_r3_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r3_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real32_r3_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real32_r3(array, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real32_r3(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r3(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real64_r3_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(3)
    real(c_double), pointer :: array_standard_bounds(:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array_standard_bounds(:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real64_r3_shape(array, shape, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_real64_r3_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r3_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real64_r3_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real64_r3(array, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real64_r3(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r3(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int32_r4_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(4)
    integer(c_int32_t), pointer :: array_standard_bounds(:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array_standard_bounds(:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int32_r4_shape(array, shape, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_int32_r4_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r4_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int32_r4_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int32_r4(array, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int32_r4(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r4(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int64_r4_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(4)
    integer(c_int64_t), pointer :: array_standard_bounds(:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array_standard_bounds(:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int64_r4_shape(array, shape, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_int64_r4_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r4_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int64_r4_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int64_r4(array, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int64_r4(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r4(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real32_r4_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(4)
    real(c_float), pointer :: array_standard_bounds(:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array_standard_bounds(:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real32_r4_shape(array, shape, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_real32_r4_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r4_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real32_r4_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real32_r4(array, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real32_r4(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r4(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real64_r4_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(4)
    real(c_double), pointer :: array_standard_bounds(:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array_standard_bounds(:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real64_r4_shape(array, shape, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_real64_r4_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r4_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real64_r4_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real64_r4(array, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real64_r4(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r4(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int32_r5_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(5)
    integer(c_int32_t), pointer :: array_standard_bounds(:,:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array_standard_bounds(:,:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int32_r5_shape(array, shape, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_int32_r5_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r5_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int32_r5_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int32_r5(array, resource)
    implicit none
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int32_r5(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r5(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_int64_r5_bounds(array, lbounds, ubounds, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(5)
    integer(c_int64_t), pointer :: array_standard_bounds(:,:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array_standard_bounds(:,:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_int64_r5_shape(array, shape, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_int64_r5_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r5_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_int64_r5_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_int64_r5(array, resource)
    implicit none
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_int64_r5(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r5(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real32_r5_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(5)
    real(c_float), pointer :: array_standard_bounds(:,:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 4
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array_standard_bounds(:,:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real32_r5_shape(array, shape, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_real32_r5_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r5_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real32_r5_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real32_r5(array, resource)
    implicit none
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real32_r5(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r5(array, resource)
    call pluto_set_label(previous_label)
end subroutine


subroutine pluto_allocate_real64_r5_bounds(array, lbounds, ubounds, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_int) :: shape(5)
    real(c_double), pointer :: array_standard_bounds(:,:,:,:,:)
    shape = ubounds - lbounds + 1
    bytes = array_size(lbounds, ubounds) * 8
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array_standard_bounds, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array_standard_bounds(:,:,:,:,:)
    else
        nullify(array)
    endif
end subroutine

subroutine pluto_allocate_real64_r5_shape(array, shape, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    integer(c_int) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_real64_r5_bounds(array, standard_lbounds, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r5_bounds(label, array, lbounds, ubounds, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_allocate_label_real64_r5_shape(label, array, shape, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end subroutine

subroutine pluto_deallocate_real64_r5(array, resource)
    implicit none
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end subroutine

subroutine pluto_deallocate_label_real64_r5(label, array, resource)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r5(array, resource)
    call pluto_set_label(previous_label)
end subroutine


end module
