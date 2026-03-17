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

use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
use pluto_module_memory_resource, only : pluto_memory_resource, pluto_get_label, pluto_set_label
implicit none
private

public :: pluto_allocate, pluto_deallocate

interface pluto_allocate
    module procedure pluto_allocate_int32_r1_bounds
    module procedure pluto_allocate_int64_r1_bounds
    module procedure pluto_allocate_real32_r1_bounds
    module procedure pluto_allocate_real64_r1_bounds
    module procedure pluto_allocate_int32_r2_bounds
    module procedure pluto_allocate_int64_r2_bounds
    module procedure pluto_allocate_real32_r2_bounds
    module procedure pluto_allocate_real64_r2_bounds
    module procedure pluto_allocate_int32_r3_bounds
    module procedure pluto_allocate_int64_r3_bounds
    module procedure pluto_allocate_real32_r3_bounds
    module procedure pluto_allocate_real64_r3_bounds
    module procedure pluto_allocate_int32_r4_bounds
    module procedure pluto_allocate_int64_r4_bounds
    module procedure pluto_allocate_real32_r4_bounds
    module procedure pluto_allocate_real64_r4_bounds
    module procedure pluto_allocate_int32_r5_bounds
    module procedure pluto_allocate_int64_r5_bounds
    module procedure pluto_allocate_real32_r5_bounds
    module procedure pluto_allocate_real64_r5_bounds
    module procedure pluto_allocate_int32_r1_shape
    module procedure pluto_allocate_int64_r1_shape
    module procedure pluto_allocate_real32_r1_shape
    module procedure pluto_allocate_real64_r1_shape
    module procedure pluto_allocate_int32_r2_shape
    module procedure pluto_allocate_int64_r2_shape
    module procedure pluto_allocate_real32_r2_shape
    module procedure pluto_allocate_real64_r2_shape
    module procedure pluto_allocate_int32_r3_shape
    module procedure pluto_allocate_int64_r3_shape
    module procedure pluto_allocate_real32_r3_shape
    module procedure pluto_allocate_real64_r3_shape
    module procedure pluto_allocate_int32_r4_shape
    module procedure pluto_allocate_int64_r4_shape
    module procedure pluto_allocate_real32_r4_shape
    module procedure pluto_allocate_real64_r4_shape
    module procedure pluto_allocate_int32_r5_shape
    module procedure pluto_allocate_int64_r5_shape
    module procedure pluto_allocate_real32_r5_shape
    module procedure pluto_allocate_real64_r5_shape
    module procedure pluto_allocate_label_int32_r1_bounds
    module procedure pluto_allocate_label_int64_r1_bounds
    module procedure pluto_allocate_label_real32_r1_bounds
    module procedure pluto_allocate_label_real64_r1_bounds
    module procedure pluto_allocate_label_int32_r2_bounds
    module procedure pluto_allocate_label_int64_r2_bounds
    module procedure pluto_allocate_label_real32_r2_bounds
    module procedure pluto_allocate_label_real64_r2_bounds
    module procedure pluto_allocate_label_int32_r3_bounds
    module procedure pluto_allocate_label_int64_r3_bounds
    module procedure pluto_allocate_label_real32_r3_bounds
    module procedure pluto_allocate_label_real64_r3_bounds
    module procedure pluto_allocate_label_int32_r4_bounds
    module procedure pluto_allocate_label_int64_r4_bounds
    module procedure pluto_allocate_label_real32_r4_bounds
    module procedure pluto_allocate_label_real64_r4_bounds
    module procedure pluto_allocate_label_int32_r5_bounds
    module procedure pluto_allocate_label_int64_r5_bounds
    module procedure pluto_allocate_label_real32_r5_bounds
    module procedure pluto_allocate_label_real64_r5_bounds
    module procedure pluto_allocate_label_int32_r1_shape
    module procedure pluto_allocate_label_int64_r1_shape
    module procedure pluto_allocate_label_real32_r1_shape
    module procedure pluto_allocate_label_real64_r1_shape
    module procedure pluto_allocate_label_int32_r2_shape
    module procedure pluto_allocate_label_int64_r2_shape
    module procedure pluto_allocate_label_real32_r2_shape
    module procedure pluto_allocate_label_real64_r2_shape
    module procedure pluto_allocate_label_int32_r3_shape
    module procedure pluto_allocate_label_int64_r3_shape
    module procedure pluto_allocate_label_real32_r3_shape
    module procedure pluto_allocate_label_real64_r3_shape
    module procedure pluto_allocate_label_int32_r4_shape
    module procedure pluto_allocate_label_int64_r4_shape
    module procedure pluto_allocate_label_real32_r4_shape
    module procedure pluto_allocate_label_real64_r4_shape
    module procedure pluto_allocate_label_int32_r5_shape
    module procedure pluto_allocate_label_int64_r5_shape
    module procedure pluto_allocate_label_real32_r5_shape
    module procedure pluto_allocate_label_real64_r5_shape
end interface

interface pluto_deallocate
    module procedure pluto_deallocate_int32_r1
    module procedure pluto_deallocate_int64_r1
    module procedure pluto_deallocate_real32_r1
    module procedure pluto_deallocate_real64_r1
    module procedure pluto_deallocate_int32_r2
    module procedure pluto_deallocate_int64_r2
    module procedure pluto_deallocate_real32_r2
    module procedure pluto_deallocate_real64_r2
    module procedure pluto_deallocate_int32_r3
    module procedure pluto_deallocate_int64_r3
    module procedure pluto_deallocate_real32_r3
    module procedure pluto_deallocate_real64_r3
    module procedure pluto_deallocate_int32_r4
    module procedure pluto_deallocate_int64_r4
    module procedure pluto_deallocate_real32_r4
    module procedure pluto_deallocate_real64_r4
    module procedure pluto_deallocate_int32_r5
    module procedure pluto_deallocate_int64_r5
    module procedure pluto_deallocate_real32_r5
    module procedure pluto_deallocate_real64_r5
    module procedure pluto_deallocate_label_int32_r1
    module procedure pluto_deallocate_label_int64_r1
    module procedure pluto_deallocate_label_real32_r1
    module procedure pluto_deallocate_label_real64_r1
    module procedure pluto_deallocate_label_int32_r2
    module procedure pluto_deallocate_label_int64_r2
    module procedure pluto_deallocate_label_real32_r2
    module procedure pluto_deallocate_label_real64_r2
    module procedure pluto_deallocate_label_int32_r3
    module procedure pluto_deallocate_label_int64_r3
    module procedure pluto_deallocate_label_real32_r3
    module procedure pluto_deallocate_label_real64_r3
    module procedure pluto_deallocate_label_int32_r4
    module procedure pluto_deallocate_label_int64_r4
    module procedure pluto_deallocate_label_real32_r4
    module procedure pluto_deallocate_label_real64_r4
    module procedure pluto_deallocate_label_int32_r5
    module procedure pluto_deallocate_label_int64_r5
    module procedure pluto_deallocate_label_real32_r5
    module procedure pluto_deallocate_label_real64_r5
end interface

interface
    module subroutine pluto_allocate_int32_r1_bounds(array, lbounds, ubounds, resource)
        integer(int32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r1_shape(array, shape, resource)
        integer(int32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r1_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r1_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int32_r1(array, resource)
        integer(int32), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int32_r1(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r1_bounds(array, lbounds, ubounds, resource)
        integer(int64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r1_shape(array, shape, resource)
        integer(int64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r1_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r1_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int64_r1(array, resource)
        integer(int64), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int64_r1(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r1_bounds(array, lbounds, ubounds, resource)
        real(real32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r1_shape(array, shape, resource)
        real(real32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r1_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r1_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real32_r1(array, resource)
        real(real32), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real32_r1(label, array, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r1_bounds(array, lbounds, ubounds, resource)
        real(real64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r1_shape(array, shape, resource)
        real(real64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r1_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: lbounds(1)
        integer(int32), intent(in) :: ubounds(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r1_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:)
        integer(int32), intent(in) :: shape(1)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real64_r1(array, resource)
        real(real64), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real64_r1(label, array, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r2_bounds(array, lbounds, ubounds, resource)
        integer(int32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r2_shape(array, shape, resource)
        integer(int32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r2_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r2_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int32_r2(array, resource)
        integer(int32), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int32_r2(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r2_bounds(array, lbounds, ubounds, resource)
        integer(int64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r2_shape(array, shape, resource)
        integer(int64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r2_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r2_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int64_r2(array, resource)
        integer(int64), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int64_r2(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r2_bounds(array, lbounds, ubounds, resource)
        real(real32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r2_shape(array, shape, resource)
        real(real32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r2_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r2_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real32_r2(array, resource)
        real(real32), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real32_r2(label, array, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r2_bounds(array, lbounds, ubounds, resource)
        real(real64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r2_shape(array, shape, resource)
        real(real64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r2_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: lbounds(2)
        integer(int32), intent(in) :: ubounds(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r2_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:)
        integer(int32), intent(in) :: shape(2)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real64_r2(array, resource)
        real(real64), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real64_r2(label, array, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r3_bounds(array, lbounds, ubounds, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r3_shape(array, shape, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r3_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r3_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int32_r3(array, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int32_r3(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r3_bounds(array, lbounds, ubounds, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r3_shape(array, shape, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r3_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r3_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int64_r3(array, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int64_r3(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r3_bounds(array, lbounds, ubounds, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r3_shape(array, shape, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r3_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r3_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real32_r3(array, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real32_r3(label, array, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r3_bounds(array, lbounds, ubounds, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r3_shape(array, shape, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r3_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: lbounds(3)
        integer(int32), intent(in) :: ubounds(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r3_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:)
        integer(int32), intent(in) :: shape(3)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real64_r3(array, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real64_r3(label, array, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r4_bounds(array, lbounds, ubounds, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r4_shape(array, shape, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r4_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r4_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int32_r4(array, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int32_r4(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r4_bounds(array, lbounds, ubounds, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r4_shape(array, shape, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r4_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r4_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int64_r4(array, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int64_r4(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r4_bounds(array, lbounds, ubounds, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r4_shape(array, shape, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r4_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r4_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real32_r4(array, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real32_r4(label, array, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r4_bounds(array, lbounds, ubounds, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r4_shape(array, shape, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r4_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: lbounds(4)
        integer(int32), intent(in) :: ubounds(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r4_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:,:)
        integer(int32), intent(in) :: shape(4)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real64_r4(array, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real64_r4(label, array, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r5_bounds(array, lbounds, ubounds, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int32_r5_shape(array, shape, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r5_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int32_r5_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int32_r5(array, resource)
        integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int32_r5(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r5_bounds(array, lbounds, ubounds, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_int64_r5_shape(array, shape, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r5_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_int64_r5_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_int64_r5(array, resource)
        integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_int64_r5(label, array, resource)
        character(len=*), intent(in) :: label
        integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r5_bounds(array, lbounds, ubounds, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real32_r5_shape(array, shape, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r5_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real32_r5_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real32_r5(array, resource)
        real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real32_r5(label, array, resource)
        character(len=*), intent(in) :: label
        real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r5_bounds(array, lbounds, ubounds, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_real64_r5_shape(array, shape, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r5_bounds(label, array, lbounds, ubounds, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: lbounds(5)
        integer(int32), intent(in) :: ubounds(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_allocate_label_real64_r5_shape(label, array, shape, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
        integer(int32), intent(in) :: shape(5)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_real64_r5(array, resource)
        real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

    module subroutine pluto_deallocate_label_real64_r5(label, array, resource)
        character(len=*), intent(in) :: label
        real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
        type(pluto_memory_resource), intent(in) :: resource
    end subroutine

end interface

end module
