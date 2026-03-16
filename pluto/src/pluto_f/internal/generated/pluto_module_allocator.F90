! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_allocator

use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
use pluto_module_memory_resource, only : pluto_memory_resource

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
    procedure, private :: pluto_allocator_allocate_int64_r1_shape
    procedure, private :: pluto_allocator_allocate_int64_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r1_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r1_bounds
    procedure, private :: pluto_allocator_allocate_real32_r1_shape
    procedure, private :: pluto_allocator_allocate_real32_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r1_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r1_bounds
    procedure, private :: pluto_allocator_allocate_real64_r1_shape
    procedure, private :: pluto_allocator_allocate_real64_r1_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r1_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r1_bounds
    procedure, private :: pluto_allocator_allocate_int32_r2_shape
    procedure, private :: pluto_allocator_allocate_int32_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r2_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r2_bounds
    procedure, private :: pluto_allocator_allocate_int64_r2_shape
    procedure, private :: pluto_allocator_allocate_int64_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r2_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r2_bounds
    procedure, private :: pluto_allocator_allocate_real32_r2_shape
    procedure, private :: pluto_allocator_allocate_real32_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r2_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r2_bounds
    procedure, private :: pluto_allocator_allocate_real64_r2_shape
    procedure, private :: pluto_allocator_allocate_real64_r2_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r2_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r2_bounds
    procedure, private :: pluto_allocator_allocate_int32_r3_shape
    procedure, private :: pluto_allocator_allocate_int32_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r3_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r3_bounds
    procedure, private :: pluto_allocator_allocate_int64_r3_shape
    procedure, private :: pluto_allocator_allocate_int64_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r3_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r3_bounds
    procedure, private :: pluto_allocator_allocate_real32_r3_shape
    procedure, private :: pluto_allocator_allocate_real32_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r3_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r3_bounds
    procedure, private :: pluto_allocator_allocate_real64_r3_shape
    procedure, private :: pluto_allocator_allocate_real64_r3_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r3_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r3_bounds
    procedure, private :: pluto_allocator_allocate_int32_r4_shape
    procedure, private :: pluto_allocator_allocate_int32_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r4_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r4_bounds
    procedure, private :: pluto_allocator_allocate_int64_r4_shape
    procedure, private :: pluto_allocator_allocate_int64_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r4_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r4_bounds
    procedure, private :: pluto_allocator_allocate_real32_r4_shape
    procedure, private :: pluto_allocator_allocate_real32_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r4_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r4_bounds
    procedure, private :: pluto_allocator_allocate_real64_r4_shape
    procedure, private :: pluto_allocator_allocate_real64_r4_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r4_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r4_bounds
    procedure, private :: pluto_allocator_allocate_int32_r5_shape
    procedure, private :: pluto_allocator_allocate_int32_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_int32_r5_shape
    procedure, private :: pluto_allocator_allocate_label_int32_r5_bounds
    procedure, private :: pluto_allocator_allocate_int64_r5_shape
    procedure, private :: pluto_allocator_allocate_int64_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_int64_r5_shape
    procedure, private :: pluto_allocator_allocate_label_int64_r5_bounds
    procedure, private :: pluto_allocator_allocate_real32_r5_shape
    procedure, private :: pluto_allocator_allocate_real32_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_real32_r5_shape
    procedure, private :: pluto_allocator_allocate_label_real32_r5_bounds
    procedure, private :: pluto_allocator_allocate_real64_r5_shape
    procedure, private :: pluto_allocator_allocate_real64_r5_bounds
    procedure, private :: pluto_allocator_allocate_label_real64_r5_shape
    procedure, private :: pluto_allocator_allocate_label_real64_r5_bounds

    generic, public :: allocate => pluto_allocator_allocate_int32_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_int32_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_int64_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_int64_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_real32_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_real32_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_real64_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_real64_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r1_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r1_bounds
    generic, public :: allocate => pluto_allocator_allocate_int32_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_int32_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_int64_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_int64_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_real32_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_real32_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_real64_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_real64_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r2_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r2_bounds
    generic, public :: allocate => pluto_allocator_allocate_int32_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_int32_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_int64_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_int64_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_real32_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_real32_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_real64_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_real64_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r3_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r3_bounds
    generic, public :: allocate => pluto_allocator_allocate_int32_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_int32_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_int64_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_int64_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_real32_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_real32_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_real64_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_real64_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r4_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r4_bounds
    generic, public :: allocate => pluto_allocator_allocate_int32_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_int32_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int32_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_int64_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_int64_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_label_int64_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_real32_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_real32_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real32_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_real64_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_real64_r5_bounds
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r5_shape
    generic, public :: allocate => pluto_allocator_allocate_label_real64_r5_bounds

    procedure, private :: pluto_allocator_deallocate_int32_r1
    procedure, private :: pluto_allocator_deallocate_label_int32_r1
    procedure, private :: pluto_allocator_deallocate_int64_r1
    procedure, private :: pluto_allocator_deallocate_label_int64_r1
    procedure, private :: pluto_allocator_deallocate_real32_r1
    procedure, private :: pluto_allocator_deallocate_label_real32_r1
    procedure, private :: pluto_allocator_deallocate_real64_r1
    procedure, private :: pluto_allocator_deallocate_label_real64_r1
    procedure, private :: pluto_allocator_deallocate_int32_r2
    procedure, private :: pluto_allocator_deallocate_label_int32_r2
    procedure, private :: pluto_allocator_deallocate_int64_r2
    procedure, private :: pluto_allocator_deallocate_label_int64_r2
    procedure, private :: pluto_allocator_deallocate_real32_r2
    procedure, private :: pluto_allocator_deallocate_label_real32_r2
    procedure, private :: pluto_allocator_deallocate_real64_r2
    procedure, private :: pluto_allocator_deallocate_label_real64_r2
    procedure, private :: pluto_allocator_deallocate_int32_r3
    procedure, private :: pluto_allocator_deallocate_label_int32_r3
    procedure, private :: pluto_allocator_deallocate_int64_r3
    procedure, private :: pluto_allocator_deallocate_label_int64_r3
    procedure, private :: pluto_allocator_deallocate_real32_r3
    procedure, private :: pluto_allocator_deallocate_label_real32_r3
    procedure, private :: pluto_allocator_deallocate_real64_r3
    procedure, private :: pluto_allocator_deallocate_label_real64_r3
    procedure, private :: pluto_allocator_deallocate_int32_r4
    procedure, private :: pluto_allocator_deallocate_label_int32_r4
    procedure, private :: pluto_allocator_deallocate_int64_r4
    procedure, private :: pluto_allocator_deallocate_label_int64_r4
    procedure, private :: pluto_allocator_deallocate_real32_r4
    procedure, private :: pluto_allocator_deallocate_label_real32_r4
    procedure, private :: pluto_allocator_deallocate_real64_r4
    procedure, private :: pluto_allocator_deallocate_label_real64_r4
    procedure, private :: pluto_allocator_deallocate_int32_r5
    procedure, private :: pluto_allocator_deallocate_label_int32_r5
    procedure, private :: pluto_allocator_deallocate_int64_r5
    procedure, private :: pluto_allocator_deallocate_label_int64_r5
    procedure, private :: pluto_allocator_deallocate_real32_r5
    procedure, private :: pluto_allocator_deallocate_label_real32_r5
    procedure, private :: pluto_allocator_deallocate_real64_r5
    procedure, private :: pluto_allocator_deallocate_label_real64_r5

    generic, public :: deallocate => pluto_allocator_deallocate_int32_r1
    generic, public :: deallocate => pluto_allocator_deallocate_label_int32_r1
    generic, public :: deallocate => pluto_allocator_deallocate_int64_r1
    generic, public :: deallocate => pluto_allocator_deallocate_label_int64_r1
    generic, public :: deallocate => pluto_allocator_deallocate_real32_r1
    generic, public :: deallocate => pluto_allocator_deallocate_label_real32_r1
    generic, public :: deallocate => pluto_allocator_deallocate_real64_r1
    generic, public :: deallocate => pluto_allocator_deallocate_label_real64_r1
    generic, public :: deallocate => pluto_allocator_deallocate_int32_r2
    generic, public :: deallocate => pluto_allocator_deallocate_label_int32_r2
    generic, public :: deallocate => pluto_allocator_deallocate_int64_r2
    generic, public :: deallocate => pluto_allocator_deallocate_label_int64_r2
    generic, public :: deallocate => pluto_allocator_deallocate_real32_r2
    generic, public :: deallocate => pluto_allocator_deallocate_label_real32_r2
    generic, public :: deallocate => pluto_allocator_deallocate_real64_r2
    generic, public :: deallocate => pluto_allocator_deallocate_label_real64_r2
    generic, public :: deallocate => pluto_allocator_deallocate_int32_r3
    generic, public :: deallocate => pluto_allocator_deallocate_label_int32_r3
    generic, public :: deallocate => pluto_allocator_deallocate_int64_r3
    generic, public :: deallocate => pluto_allocator_deallocate_label_int64_r3
    generic, public :: deallocate => pluto_allocator_deallocate_real32_r3
    generic, public :: deallocate => pluto_allocator_deallocate_label_real32_r3
    generic, public :: deallocate => pluto_allocator_deallocate_real64_r3
    generic, public :: deallocate => pluto_allocator_deallocate_label_real64_r3
    generic, public :: deallocate => pluto_allocator_deallocate_int32_r4
    generic, public :: deallocate => pluto_allocator_deallocate_label_int32_r4
    generic, public :: deallocate => pluto_allocator_deallocate_int64_r4
    generic, public :: deallocate => pluto_allocator_deallocate_label_int64_r4
    generic, public :: deallocate => pluto_allocator_deallocate_real32_r4
    generic, public :: deallocate => pluto_allocator_deallocate_label_real32_r4
    generic, public :: deallocate => pluto_allocator_deallocate_real64_r4
    generic, public :: deallocate => pluto_allocator_deallocate_label_real64_r4
    generic, public :: deallocate => pluto_allocator_deallocate_int32_r5
    generic, public :: deallocate => pluto_allocator_deallocate_label_int32_r5
    generic, public :: deallocate => pluto_allocator_deallocate_int64_r5
    generic, public :: deallocate => pluto_allocator_deallocate_label_int64_r5
    generic, public :: deallocate => pluto_allocator_deallocate_real32_r5
    generic, public :: deallocate => pluto_allocator_deallocate_label_real32_r5
    generic, public :: deallocate => pluto_allocator_deallocate_real64_r5
    generic, public :: deallocate => pluto_allocator_deallocate_label_real64_r5
end type

interface pluto_make_allocator
    module procedure pluto_make_allocator_type
    module procedure pluto_make_allocator_name
end interface

interface
module function pluto_make_allocator_type(resource) result(allocator)
    type(pluto_allocator) :: allocator
    type(pluto_memory_resource) :: resource
end function

module function pluto_make_allocator_name(resource) result(allocator)
    type(pluto_allocator) :: allocator
    character(len=*), target, intent(in) :: resource
end function
end interface

interface
module subroutine pluto_allocator_allocate_int32_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_int32_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_int64_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_int64_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_real32_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_real32_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_real64_r1_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_real64_r1_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r1_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: lbounds(1), ubounds(1)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r1_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:)
    integer(int32), intent(in) :: shape(1)
end subroutine
module subroutine pluto_allocator_allocate_int32_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_int32_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_int64_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_int64_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_real32_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_real32_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_real64_r2_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_real64_r2_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r2_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: lbounds(2), ubounds(2)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r2_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:)
    integer(int32), intent(in) :: shape(2)
end subroutine
module subroutine pluto_allocator_allocate_int32_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_int32_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_int64_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_int64_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_real32_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_real32_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_real64_r3_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_real64_r3_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r3_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: lbounds(3), ubounds(3)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r3_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:)
    integer(int32), intent(in) :: shape(3)
end subroutine
module subroutine pluto_allocator_allocate_int32_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_int32_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_int64_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_int64_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_real32_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_real32_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_real64_r4_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_real64_r4_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r4_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: lbounds(4), ubounds(4)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r4_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:,:)
    integer(int32), intent(in) :: shape(4)
end subroutine
module subroutine pluto_allocator_allocate_int32_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_int32_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_label_int32_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_int64_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_int64_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_label_int64_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_real32_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_real32_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_label_real32_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_real64_r5_bounds(this, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_real64_r5_shape(this, array, shape)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r5_bounds(this, label, array, lbounds, ubounds)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: lbounds(5), ubounds(5)
end subroutine
module subroutine pluto_allocator_allocate_label_real64_r5_shape(this, label, array, shape)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(int32), intent(in) :: shape(5)
end subroutine

module subroutine pluto_allocator_deallocate_int32_r1(this, array)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int32_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_int64_r1(this, array)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int64_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_real32_r1(this, array)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real32_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_real64_r1(this, array)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real64_r1(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:)
end subroutine
module subroutine pluto_allocator_deallocate_int32_r2(this, array)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int32_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int64_r2(this, array)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int64_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real32_r2(this, array)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real32_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real64_r2(this, array)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real64_r2(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int32_r3(this, array)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int32_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int64_r3(this, array)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int64_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real32_r3(this, array)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real32_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real64_r3(this, array)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real64_r3(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int32_r4(this, array)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int32_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int64_r4(this, array)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int64_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real32_r4(this, array)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real32_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real64_r4(this, array)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real64_r4(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int32_r5(this, array)
    class(pluto_allocator) :: this
    integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int32_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int32), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_int64_r5(this, array)
    class(pluto_allocator) :: this
    integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_int64_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    integer(int64), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real32_r5(this, array)
    class(pluto_allocator) :: this
    real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real32_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real32), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_real64_r5(this, array)
    class(pluto_allocator) :: this
    real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine
module subroutine pluto_allocator_deallocate_label_real64_r5(this, label, array)
    class(pluto_allocator) :: this
    character(len=*), intent(in) :: label
    real(real64), pointer, intent(inout) :: array(:,:,:,:,:)
end subroutine

end interface
end module
