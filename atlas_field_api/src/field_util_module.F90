! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module field_util_module

! WARNING: This module is not part of the public API!

private
public :: field_threads, field_my_thread
public :: optional_arg
public :: set_host_resource
public :: set_device_mapped_host_resource

contains

function field_threads()
#ifdef _OPENMP
    use omp_lib
    integer :: field_threads
    field_threads = omp_get_max_threads()
#else
    integer :: field_threads
    field_threads = 1
#endif
end function

function field_my_thread()
#ifdef _OPENMP
    use omp_lib
    integer :: field_my_thread
    field_my_thread = omp_get_thread_num()
#else
    integer :: field_my_thread
    field_my_thread = 1
#endif
end function

function field_pooled_default()
    logical :: field_pooled_default
    field_pooled_default = .false.
end function

function field_pinned_default()
    logical :: field_pinned_default
    field_pinned_default = .false.
end function

function optional_arg(VALUE, DEFAULT)
    logical :: optional_arg
    logical, intent(in), optional :: VALUE
    logical :: DEFAULT
    if (present(VALUE)) then
      optional_arg = VALUE
    else
      optional_arg = DEFAULT
    endif
end function

subroutine set_device_mapped_host_resource()
    ! Method to select which resource to use for non-contiguous host-wrapped fields for the device_mapped_host_data pointer with atlas/array/native/NativeDataStore.h
    ! If this is not set, then pluto%host_resource() is used by default
    ! We could play with this and if not OK we can use the buddy_allocator as in Field-API
    use pluto_module
    if (.not. pluto%has_registered_resource("atlas::array::device_mapped_host_resource")) then
        call pluto%register_resource("atlas::array::device_mapped_host_resource", pluto%host_pool_resource())
    endif
end subroutine

subroutine set_host_resource(PINNED, POOLED)
    ! Method to select which host resource to use, based on values of PINNED and POOLED
    use pluto_module
    logical, intent(in), optional :: PINNED
    logical, intent(in), optional :: POOLED
    logical :: _pooled
    logical :: _pinned
    _pooled = optional_arg(POOLED, DEFAULT=field_pooled_default())
    _pinned = optional_arg(PINNED, DEfAULT=field_pinned_default())
    if (_pooled .and. _pinned) then
        call pluto%host%set_default_resource(pluto%pinned_pool_resource())
    else if (_pooled .and. .not. _pinned) then
        call pluto%host%set_default_resource(pluto%host_pool_resource())
    else if (.not. _pooled .and. _pinned) then
        call pluto%host%set_default_resource(pluto%pinned_resource())
    else if (.not. _pooled .and. .not. _pinned) then
        call pluto%host%set_default_resource(pluto%host_resource())
    endif
end subroutine

end module
