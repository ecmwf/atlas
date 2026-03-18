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
                                              & pluto_register_memory_resource_adaptor, &
                                              & pluto_reserve
use pluto_module_allocator,              only : pluto_allocator
use pluto_module_host,                   only : pluto_host_t
use pluto_module_device,                 only : pluto_device_t
use pluto_module_runtime,                only : pluto_devices
use pluto_module_scope,                  only : pluto_scope_t
use pluto_module_trace,                  only : pluto_trace_t
use pluto_module_memory,                 only : pluto_memory_t
use pluto_module_mpi,                    only : pluto_mpi_t
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

type pluto_t
    type(pluto_host_t)   :: host
    type(pluto_device_t) :: device
    type(pluto_scope_t)  :: scope
    type(pluto_trace_t)  :: trace
    type(pluto_mpi_t)    :: mpi
    type(pluto_memory_t) :: memory
contains
    procedure, nopass :: devices => pluto_devices

    procedure, nopass :: reserve => pluto_reserve

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

    procedure, private, nopass :: pluto_release_all
    procedure, private, nopass :: pluto_release_resource
    generic :: release => pluto_release_all, pluto_release_resource

end type

type(pluto_t) :: pluto

contains

function make_allocator_type(resource) result(allocator)
    type(pluto_allocator) :: allocator
    type(pluto_memory_resource), intent(in) :: resource
    call allocator%init(resource)
end function
    
function make_allocator_name(resource) result(allocator)
    type(pluto_allocator) :: allocator
    character(len=*), intent(in) :: resource
    call allocator%init(resource)
end function

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

end module
