! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

submodule(pluto_module_allocator) pluto_module_allocator_smod

use, intrinsic :: iso_c_binding      , only : c_associated
use pluto_module_abort               , only : pluto_abort
use pluto_module_memory_resource     , only : pluto_get_registered_resource
use pluto_module_allocate_deallocate , only : pluto_allocate, pluto_deallocate

implicit none


contains

subroutine assert_allocator_is_setup(allocator)
    type(pluto_allocator) :: allocator
    if (.not. c_associated(allocator%memory_resource%c_memory_resource)) then
        call pluto_abort("pluto_allocator has not been assigned or setup properly.&
          & Please ensure that the allocator is created by pluto%make_allocator or&
          & pluto%{host,device}%make_allocator.")
    endif
end subroutine

module procedure pluto_make_allocator_type
    allocator%memory_resource%c_memory_resource = resource%c_memory_resource
end procedure

module procedure pluto_make_allocator_name
    allocator%memory_resource = pluto_get_registered_resource(resource)
end procedure

module procedure pluto_allocator_allocate_int32_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int32_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int32_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int64_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int64_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real32_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real32_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r1_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r1_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real64_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real64_r1
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int32_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int32_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int64_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int64_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real32_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real32_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r2_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r2_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real64_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real64_r2
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int32_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int32_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int64_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int64_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real32_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real32_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r3_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r3_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real64_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real64_r3
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int32_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int32_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int64_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int64_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real32_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real32_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r4_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r4_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real64_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real64_r4
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int32_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int32_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int32_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int32_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_int64_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_int64_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_int64_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_int64_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real32_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real32_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real32_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real32_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_real64_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r5_bounds
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, lbounds, ubounds, this%memory_resource)
end procedure

module procedure pluto_allocator_allocate_label_real64_r5_shape
    call assert_allocator_is_setup(this)
    call pluto_allocate(label, array, shape, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_real64_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(array, this%memory_resource)
end procedure

module procedure pluto_allocator_deallocate_label_real64_r5
    call assert_allocator_is_setup(this)
    call pluto_deallocate(label, array, this%memory_resource)
end procedure


end submodule
