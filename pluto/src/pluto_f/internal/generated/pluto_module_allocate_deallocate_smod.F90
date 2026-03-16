! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

submodule(pluto_module_allocate_deallocate) pluto_module_allocate_deallocate

use, intrinsic :: iso_c_binding, only : c_loc, c_ptr, c_size_t, c_f_pointer
implicit none


contains


module procedure pluto_allocate_int32_r1_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(1)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int32_r1_shape
    integer(int32) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_int32_r1_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int32_r1_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int32_r1_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int32_r1
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int32_r1
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r1(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int64_r1_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(1)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int64_r1_shape
    integer(int32) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_int64_r1_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int64_r1_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int64_r1_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int64_r1
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int64_r1
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r1(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real32_r1_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(1)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real32_r1_shape
    integer(int32) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_real32_r1_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real32_r1_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real32_r1_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real32_r1
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real32_r1
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r1(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real64_r1_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(1)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real64_r1_shape
    integer(int32) :: standard_lbounds(1)
    standard_lbounds = 1
    call pluto_allocate_real64_r1_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real64_r1_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r1_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real64_r1_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r1_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real64_r1
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real64_r1
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r1(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int32_r2_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(2)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int32_r2_shape
    integer(int32) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_int32_r2_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int32_r2_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int32_r2_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int32_r2
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int32_r2
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r2(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int64_r2_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(2)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int64_r2_shape
    integer(int32) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_int64_r2_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int64_r2_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int64_r2_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int64_r2
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int64_r2
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r2(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real32_r2_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(2)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real32_r2_shape
    integer(int32) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_real32_r2_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real32_r2_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real32_r2_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real32_r2
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real32_r2
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r2(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real64_r2_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(2)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real64_r2_shape
    integer(int32) :: standard_lbounds(2)
    standard_lbounds = 1
    call pluto_allocate_real64_r2_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real64_r2_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r2_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real64_r2_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r2_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real64_r2
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real64_r2
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r2(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int32_r3_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(3)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int32_r3_shape
    integer(int32) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_int32_r3_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int32_r3_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int32_r3_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int32_r3
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int32_r3
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r3(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int64_r3_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(3)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int64_r3_shape
    integer(int32) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_int64_r3_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int64_r3_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int64_r3_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int64_r3
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int64_r3
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r3(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real32_r3_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(3)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real32_r3_shape
    integer(int32) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_real32_r3_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real32_r3_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real32_r3_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real32_r3
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real32_r3
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r3(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real64_r3_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(3)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real64_r3_shape
    integer(int32) :: standard_lbounds(3)
    standard_lbounds = 1
    call pluto_allocate_real64_r3_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real64_r3_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r3_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real64_r3_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r3_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real64_r3
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real64_r3
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r3(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int32_r4_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(4)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int32_r4_shape
    integer(int32) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_int32_r4_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int32_r4_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int32_r4_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int32_r4
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int32_r4
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r4(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int64_r4_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(4)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int64_r4_shape
    integer(int32) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_int64_r4_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int64_r4_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int64_r4_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int64_r4
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int64_r4
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r4(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real32_r4_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(4)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real32_r4_shape
    integer(int32) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_real32_r4_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real32_r4_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real32_r4_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real32_r4
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real32_r4
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r4(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real64_r4_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(4)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real64_r4_shape
    integer(int32) :: standard_lbounds(4)
    standard_lbounds = 1
    call pluto_allocate_real64_r4_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real64_r4_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r4_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real64_r4_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r4_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real64_r4
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real64_r4
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r4(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int32_r5_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(5)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
        bytes = bytes * shape(5)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int32_r5_shape
    integer(int32) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_int32_r5_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int32_r5_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int32_r5_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int32_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int32_r5
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int32_r5
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int32_r5(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_int64_r5_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(5)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
        bytes = bytes * shape(5)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_int64_r5_shape
    integer(int32) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_int64_r5_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_int64_r5_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_int64_r5_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_int64_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_int64_r5
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_int64_r5
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_int64_r5(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real32_r5_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(5)
    shape = ubounds - lbounds + 1
    bytes = 4
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
        bytes = bytes * shape(5)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real32_r5_shape
    integer(int32) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_real32_r5_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real32_r5_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real32_r5_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real32_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real32_r5
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 4
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real32_r5
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real32_r5(array, resource)
    call pluto_set_label(previous_label)
end procedure


module procedure pluto_allocate_real64_r5_bounds
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    integer(c_size_t) :: shape(5)
    shape = ubounds - lbounds + 1
    bytes = 8
        bytes = bytes * shape(1)
        bytes = bytes * shape(2)
        bytes = bytes * shape(3)
        bytes = bytes * shape(4)
        bytes = bytes * shape(5)
    if (bytes > 0) then
        call resource%allocate(mem, bytes)
        call c_f_pointer(mem, array, shape)
        array(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => array
    else
        nullify(array)
    endif
end procedure

module procedure pluto_allocate_real64_r5_shape
    integer(int32) :: standard_lbounds(5)
    standard_lbounds = 1
    call pluto_allocate_real64_r5_bounds(array, standard_lbounds, shape, resource)
end procedure

module procedure pluto_allocate_label_real64_r5_bounds
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r5_bounds(array, lbounds, ubounds, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_allocate_label_real64_r5_shape
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_allocate_real64_r5_shape(array, shape, resource)
    call pluto_set_label(previous_label)
end procedure

module procedure pluto_deallocate_real64_r5
    implicit none
    type(c_ptr) :: mem
    integer(c_size_t) :: bytes
    bytes = size(array) * 8
    if (bytes > 0) then
        mem = c_loc(array(lbound(array,1), lbound(array,2), lbound(array,3), lbound(array,4), lbound(array,5)))
        call resource%deallocate(mem, bytes)
    endif
    nullify(array)
end procedure

module procedure pluto_deallocate_label_real64_r5
    implicit none
    character(len=:), allocatable :: previous_label
    previous_label = pluto_get_label()
    call pluto_set_label(label)
    call pluto_deallocate_real64_r5(array, resource)
    call pluto_set_label(previous_label)
end procedure


end submodule
