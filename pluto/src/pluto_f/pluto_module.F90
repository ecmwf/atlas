! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module

use pluto_module_memory_resource,        only : pluto_memory_resource, &
                                              & pluto_has_registered_resource, &
                                              & pluto_get_registered_resource, &
                                              & pluto_register_resource, &
                                              & pluto_unregister_resource, &
                                              & pluto_new_delete_resource, &
                                              & pluto_null_memory_resource, &
                                              & pluto_host_resource, &
                                              & pluto_pinned_resource, &
                                              & pluto_device_resource, &
                                              & pluto_managed_resource, &
                                              & pluto_mpi_resource, &
                                              & pluto_host_pool_resource, &
                                              & pluto_pinned_pool_resource, &
                                              & pluto_device_pool_resource, &
                                              & pluto_managed_pool_resource, &
                                              & pluto_mpi_pool_resource, &
                                              & pluto_set_label, &
                                              & pluto_unset_label, &
                                              & pluto_get_label, &
                                              & pluto_register_memory_resource_adaptor
use pluto_module_allocator,              only : pluto_allocator
use pluto_module_host,                   only : pluto_host_t
use pluto_module_device,                 only : pluto_device_t
use pluto_module_scope,                  only : pluto_scope_t
use pluto_module_trace,                  only : pluto_trace_t
use pluto_module_allocate_deallocate,    only : pluto_allocate, pluto_deallocate

implicit none
private

public :: pluto_memory_resource
public :: pluto_allocator
public :: pluto_allocate
public :: pluto_deallocate
public :: pluto_set_label
public :: pluto_unset_label
public :: pluto_get_label
public :: pluto

type pluto_mpi_t
contains
    procedure, nopass :: init => pluto_mpi_init
    procedure, nopass :: finalize => pluto_mpi_finalize
end type

type pluto_memory_t
contains
    procedure, nopass :: report => pluto_memory_report
end type

type pluto_t
    type(pluto_host_t)   :: host
    type(pluto_device_t) :: device
    type(pluto_scope_t)  :: scope
    type(pluto_trace_t)  :: trace
    type(pluto_mpi_t)    :: mpi
    type(pluto_memory_t) :: memory
contains
    procedure, nopass :: devices => pluto_devices

    procedure, nopass :: has_registered_resource => pluto_has_registered_resource
    procedure, nopass :: get_registered_resource => pluto_get_registered_resource

    procedure, nopass :: register_resource   => pluto_register_resource
    procedure, nopass :: unregister_resource => pluto_unregister_resource

    procedure, nopass :: register_memory_resource_adaptor => pluto_register_memory_resource_adaptor

    procedure, nopass :: new_delete_resource   => pluto_new_delete_resource
    procedure, nopass :: null_memory_resource  => pluto_null_memory_resource
    procedure, nopass :: host_resource         => pluto_host_resource
    procedure, nopass :: pinned_resource       => pluto_pinned_resource
    procedure, nopass :: device_resource       => pluto_device_resource
    procedure, nopass :: managed_resource      => pluto_managed_resource
    procedure, nopass :: mpi_resource          => pluto_mpi_resource
    procedure, nopass :: host_pool_resource    => pluto_host_pool_resource
    procedure, nopass :: pinned_pool_resource  => pluto_pinned_pool_resource
    procedure, nopass :: device_pool_resource  => pluto_device_pool_resource
    procedure, nopass :: managed_pool_resource => pluto_managed_pool_resource
    procedure, nopass :: mpi_pool_resource     => pluto_mpi_pool_resource

    procedure, nopass :: set_label   => pluto_set_label
    procedure, nopass :: unset_label => pluto_unset_label
    procedure, nopass :: get_label   => pluto_get_label

    procedure, nopass, private :: make_allocator_type
    procedure, nopass, private :: make_allocator_name
    generic :: make_allocator => make_allocator_type, make_allocator_name

    procedure, nopass :: reserve => pluto_reserve

    procedure, private, nopass :: pluto_release_all
    procedure, private, nopass :: pluto_release_resource
    generic :: release => pluto_release_all, pluto_release_resource

end type

type(pluto_t) :: pluto

contains

function pluto_devices()
    use, intrinsic :: iso_fortran_env, only: int32
    integer(int32) :: pluto_devices
    interface
        function c_pluto_devices() result(devices) bind(c)
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int) :: devices
        end function
    end interface
    pluto_devices = c_pluto_devices()
end function

function make_allocator_type(resource) result(allocator)
    use pluto_module_allocator, only : pluto_allocator, pluto_make_allocator
    ! class(pluto_t), intent(in) :: this
    type(pluto_allocator) :: allocator
    type(pluto_memory_resource), intent(in) :: resource
    allocator = pluto_make_allocator(resource)
    allocator%memory_resource = resource
end function
    
function make_allocator_name(resource) result(allocator)
    use pluto_module_allocator, only : pluto_allocator, pluto_make_allocator
    ! class(pluto_t), intent(in) :: this
    type(pluto_allocator) :: allocator
    character(len=*), intent(in) :: resource
    allocator = pluto_make_allocator(resource)
end function

subroutine pluto_mpi_init()
    interface
        subroutine c_pluto_mpi_init() bind(c)
        end subroutine
    end interface
    call c_pluto_mpi_init()
end subroutine

subroutine pluto_mpi_finalize()
    interface
        subroutine c_pluto_mpi_finalize() bind(c)
        end subroutine
    end interface
    call c_pluto_mpi_finalize()
end subroutine

subroutine pluto_reserve(resource, size)
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: iso_c_binding, only: c_size_t
    implicit none
    type(pluto_memory_resource), intent(in) :: resource
    integer(int64), intent(in) :: size
    call resource%reserve(int(size,c_size_t))
end subroutine

subroutine pluto_release_all()
    interface
        subroutine c_pluto_release() bind(c)
        end subroutine
    end interface
    call c_pluto_release()
end subroutine

subroutine pluto_release_resource(resource)
    type(pluto_memory_resource), intent(in) :: resource
    call resource%release()
end subroutine

function c_ptr_to_string(str_c_ptr,str_size) result(string)
  use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_size_t, c_f_pointer
  type(c_ptr), intent(in) :: str_c_ptr
  integer(c_size_t), intent(in) :: str_size
  character(kind=c_char,len=:), allocatable :: string
  character(kind=c_char,len=1), pointer  :: str_f_ptr(:)
  integer :: c
  call c_f_pointer( str_c_ptr , str_f_ptr, [str_size] )
  allocate( character(len=(str_size)) :: string )
  do c=1,str_size
    string(c:c) = str_f_ptr(c)
  enddo
end function

function pluto_memory_report() result(string)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    character(len=:), allocatable :: string
    interface
        subroutine c_pluto_memory_report(str_c_ptr, str_size) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), intent(out) :: str_c_ptr
            integer(c_size_t), intent(out) :: str_size
        end subroutine
        subroutine c_pluto_str_delete(str_c_ptr, str_size) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value, intent(in) :: str_c_ptr
            integer(c_size_t), value, intent(in) :: str_size
        end subroutine
    end interface
    type(c_ptr) :: str_c_ptr
    integer(c_size_t) :: str_size
    call c_pluto_memory_report(str_c_ptr, str_size)
    string = c_ptr_to_string(str_c_ptr, str_size)
    call c_pluto_str_delete(str_c_ptr, str_size)
end function

end module
