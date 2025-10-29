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
public :: set_environment_variable
public :: set_device_memory
public :: get_device_memory


interface
function setenv(varname, value, flag) bind(c) result (r)
use, intrinsic :: iso_c_binding, only : c_char, c_int
character(kind=c_char), dimension(*) :: varname
character(kind=c_char), dimension(*) :: value
integer(kind=c_int), value :: flag
end function setenv
end interface

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
    use pluto_module, only : pluto
    if (.not. pluto%has_registered_resource("atlas::array::device_mapped_host_resource")) then
        call pluto%register_resource("atlas::array::device_mapped_host_resource", pluto%host_pool_resource())
    endif
end subroutine

subroutine set_host_resource(PINNED, POOLED)
    ! Method to select which host resource to use, based on values of PINNED and POOLED
    use pluto_module
    logical, intent(in), optional :: PINNED
    logical, intent(in), optional :: POOLED
    logical :: pooled_
    logical :: pinned_
    pooled_ = optional_arg(POOLED, DEFAULT=field_pooled_default())
    pinned_ = optional_arg(PINNED, DEfAULT=field_pinned_default())
    if (pooled_ .and. pinned_) then
        call pluto%host%set_default_resource(pluto%pinned_pool_resource())
    else if (pooled_ .and. .not. pinned_) then
        call pluto%host%set_default_resource(pluto%host_pool_resource())
    else if (.not. pooled_ .and. pinned_) then
        call pluto%host%set_default_resource(pluto%pinned_resource())
    else if (.not. pooled_ .and. .not. pinned_) then
        call pluto%host%set_default_resource(pluto%host_resource())
    endif
end subroutine

subroutine set_environment_variable(var_name, value, overwrite, return_code)
    use iso_c_binding, only: c_int, c_null_char
    implicit none
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in) :: value
    logical, intent(in), optional :: overwrite
    integer, intent(out), optional :: return_code

    integer(kind=c_int) :: rc, return_value, c_overwrite

    c_overwrite = 1
    if (present(overwrite)) then
        if (overwrite) then
            c_overwrite = 1
        else
            c_overwrite = 0
        endif
    endif
    return_value = setenv(trim(var_name)//c_null_char, &
    trim(value)//c_null_char, c_overwrite)

    if (present(return_code)) return_code = return_value
end subroutine set_environment_variable

function get_device_memory()
    use iso_c_binding, only : c_int
    use pluto_module, only : pluto
    logical :: get_device_memory
    get_device_memory = .false.
    if (pluto%devices() > 0) then
      get_device_memory = .true.
      return
    endif
    block
        character(len=5) :: value
        integer :: status, length
        call get_environment_variable("ATLAS_SIMULATE_DEVICE_MEMORY", value, length, status)
        if (status == 0) THEN
            if (trim(value) == "1") then
                get_device_memory=.true.
            endif
        endif
    end block
end function

subroutine set_device_memory()
    call set_environment_variable("ATLAS_SIMULATE_DEVICE_MEMORY", "1")
end subroutine

end module
