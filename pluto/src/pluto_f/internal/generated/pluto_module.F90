! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module

use, intrinsic :: iso_c_binding,         only : c_int, c_int32_t, c_int64_t,c_double, c_float
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
                                              & pluto_memory_pool_resource_reserve_int32, &
                                              & pluto_memory_pool_resource_reserve_int64, &
                                              & pluto_memory_pool_resource_reserve_real32, &
                                              & pluto_memory_pool_resource_reserve_real64, &
                                              & pluto_memory_pool_resource_release, &
                                              & pluto_register_memory_resource_adaptor
use pluto_module_allocator,              only : pluto_allocator
use pluto_module_host,                   only : pluto_host_t
use pluto_module_device,                 only : pluto_device_t
use pluto_module_scope,                  only : pluto_scope_t
use pluto_module_trace,                  only : pluto_trace_t
use pluto_module_allocate_deallocate,    only : pluto_allocate, pluto_deallocate

implicit none
private


public :: pluto, pluto_memory_resource, pluto_allocator

type pluto_mpi_t
contains
    procedure, nopass :: init => pluto_mpi_init
    procedure, nopass :: finalize => pluto_mpi_finalize
end type

type pluto_t
    type(pluto_host_t)   :: host
    type(pluto_device_t) :: device
    type(pluto_scope_t)  :: scope
    type(pluto_trace_t)  :: trace
    type(pluto_mpi_t)    :: mpi
contains
    procedure, nopass :: devices => pluto_devices
    procedure, nopass :: has_registered_resource => pluto_has_registered_resource
    procedure, nopass :: get_registered_resource => pluto_get_registered_resource
    procedure, nopass :: register_resource => pluto_register_resource
    procedure, nopass :: unregister_resource => pluto_unregister_resource
    procedure, nopass :: register_memory_resource_adaptor => pluto_register_memory_resource_adaptor
    procedure, nopass :: new_delete_resource => pluto_new_delete_resource
    procedure, nopass :: null_memory_resource => pluto_null_memory_resource
    procedure, nopass :: host_resource => pluto_host_resource
    procedure, nopass :: pinned_resource => pluto_pinned_resource
    procedure, nopass :: device_resource => pluto_device_resource
    procedure, nopass :: managed_resource => pluto_managed_resource
    procedure, nopass :: mpi_resource => pluto_mpi_resource
    procedure, nopass :: host_pool_resource => pluto_host_pool_resource
    procedure, nopass :: pinned_pool_resource => pluto_pinned_pool_resource
    procedure, nopass :: device_pool_resource => pluto_device_pool_resource
    procedure, nopass :: managed_pool_resource => pluto_managed_pool_resource
    procedure, nopass :: mpi_pool_resource => pluto_mpi_pool_resource

    procedure, nopass :: set_label   => pluto_set_label
    procedure, nopass :: unset_label => pluto_unset_label
    procedure, nopass :: get_label   => pluto_get_label

    procedure, nopass, private :: make_allocator_type
    procedure, nopass, private :: make_allocator_name
    generic :: make_allocator => make_allocator_type, make_allocator_name
    
    procedure, private ,nopass :: reserve_int32 => pluto_memory_pool_resource_reserve_int32
    generic :: reserve => reserve_int32
    procedure, private ,nopass :: reserve_int64 => pluto_memory_pool_resource_reserve_int64
    generic :: reserve => reserve_int64
    procedure, private ,nopass :: reserve_real32 => pluto_memory_pool_resource_reserve_real32
    generic :: reserve => reserve_real32
    procedure, private ,nopass :: reserve_real64 => pluto_memory_pool_resource_reserve_real64
    generic :: reserve => reserve_real64
    
    procedure, nopass :: release => pluto_memory_pool_resource_release

    procedure, nopass, private :: pluto_allocate_int32_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_int32_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_int32_r1_shape_type
    procedure, nopass, private :: pluto_allocate_int32_r1_shape_name
    procedure, nopass, private :: pluto_allocate_label_int32_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int32_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int32_r1_shape_type
    procedure, nopass, private :: pluto_allocate_label_int32_r1_shape_name
    generic :: allocate => &
        & pluto_allocate_int32_r1_bounds_type, &
        & pluto_allocate_int32_r1_bounds_name, &
        & pluto_allocate_int32_r1_shape_type, &
        & pluto_allocate_int32_r1_shape_name, &
        & pluto_allocate_label_int32_r1_bounds_type, &
        & pluto_allocate_label_int32_r1_bounds_name, &
        & pluto_allocate_label_int32_r1_shape_type, &
        & pluto_allocate_label_int32_r1_shape_name
    procedure, nopass, private :: pluto_deallocate_int32_r1_type
    procedure, nopass, private :: pluto_deallocate_int32_r1_name
    procedure, nopass, private :: pluto_deallocate_label_int32_r1_type
    procedure, nopass, private :: pluto_deallocate_label_int32_r1_name
    generic :: deallocate => &
        & pluto_deallocate_int32_r1_type, &
        & pluto_deallocate_int32_r1_name, &
        & pluto_deallocate_label_int32_r1_type, &
        & pluto_deallocate_label_int32_r1_name
    procedure, nopass, private :: pluto_allocate_int64_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_int64_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_int64_r1_shape_type
    procedure, nopass, private :: pluto_allocate_int64_r1_shape_name
    procedure, nopass, private :: pluto_allocate_label_int64_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int64_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int64_r1_shape_type
    procedure, nopass, private :: pluto_allocate_label_int64_r1_shape_name
    generic :: allocate => &
        & pluto_allocate_int64_r1_bounds_type, &
        & pluto_allocate_int64_r1_bounds_name, &
        & pluto_allocate_int64_r1_shape_type, &
        & pluto_allocate_int64_r1_shape_name, &
        & pluto_allocate_label_int64_r1_bounds_type, &
        & pluto_allocate_label_int64_r1_bounds_name, &
        & pluto_allocate_label_int64_r1_shape_type, &
        & pluto_allocate_label_int64_r1_shape_name
    procedure, nopass, private :: pluto_deallocate_int64_r1_type
    procedure, nopass, private :: pluto_deallocate_int64_r1_name
    procedure, nopass, private :: pluto_deallocate_label_int64_r1_type
    procedure, nopass, private :: pluto_deallocate_label_int64_r1_name
    generic :: deallocate => &
        & pluto_deallocate_int64_r1_type, &
        & pluto_deallocate_int64_r1_name, &
        & pluto_deallocate_label_int64_r1_type, &
        & pluto_deallocate_label_int64_r1_name
    procedure, nopass, private :: pluto_allocate_real32_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_real32_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_real32_r1_shape_type
    procedure, nopass, private :: pluto_allocate_real32_r1_shape_name
    procedure, nopass, private :: pluto_allocate_label_real32_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real32_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real32_r1_shape_type
    procedure, nopass, private :: pluto_allocate_label_real32_r1_shape_name
    generic :: allocate => &
        & pluto_allocate_real32_r1_bounds_type, &
        & pluto_allocate_real32_r1_bounds_name, &
        & pluto_allocate_real32_r1_shape_type, &
        & pluto_allocate_real32_r1_shape_name, &
        & pluto_allocate_label_real32_r1_bounds_type, &
        & pluto_allocate_label_real32_r1_bounds_name, &
        & pluto_allocate_label_real32_r1_shape_type, &
        & pluto_allocate_label_real32_r1_shape_name
    procedure, nopass, private :: pluto_deallocate_real32_r1_type
    procedure, nopass, private :: pluto_deallocate_real32_r1_name
    procedure, nopass, private :: pluto_deallocate_label_real32_r1_type
    procedure, nopass, private :: pluto_deallocate_label_real32_r1_name
    generic :: deallocate => &
        & pluto_deallocate_real32_r1_type, &
        & pluto_deallocate_real32_r1_name, &
        & pluto_deallocate_label_real32_r1_type, &
        & pluto_deallocate_label_real32_r1_name
    procedure, nopass, private :: pluto_allocate_real64_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_real64_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_real64_r1_shape_type
    procedure, nopass, private :: pluto_allocate_real64_r1_shape_name
    procedure, nopass, private :: pluto_allocate_label_real64_r1_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real64_r1_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real64_r1_shape_type
    procedure, nopass, private :: pluto_allocate_label_real64_r1_shape_name
    generic :: allocate => &
        & pluto_allocate_real64_r1_bounds_type, &
        & pluto_allocate_real64_r1_bounds_name, &
        & pluto_allocate_real64_r1_shape_type, &
        & pluto_allocate_real64_r1_shape_name, &
        & pluto_allocate_label_real64_r1_bounds_type, &
        & pluto_allocate_label_real64_r1_bounds_name, &
        & pluto_allocate_label_real64_r1_shape_type, &
        & pluto_allocate_label_real64_r1_shape_name
    procedure, nopass, private :: pluto_deallocate_real64_r1_type
    procedure, nopass, private :: pluto_deallocate_real64_r1_name
    procedure, nopass, private :: pluto_deallocate_label_real64_r1_type
    procedure, nopass, private :: pluto_deallocate_label_real64_r1_name
    generic :: deallocate => &
        & pluto_deallocate_real64_r1_type, &
        & pluto_deallocate_real64_r1_name, &
        & pluto_deallocate_label_real64_r1_type, &
        & pluto_deallocate_label_real64_r1_name
    procedure, nopass, private :: pluto_allocate_int32_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_int32_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_int32_r2_shape_type
    procedure, nopass, private :: pluto_allocate_int32_r2_shape_name
    procedure, nopass, private :: pluto_allocate_label_int32_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int32_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int32_r2_shape_type
    procedure, nopass, private :: pluto_allocate_label_int32_r2_shape_name
    generic :: allocate => &
        & pluto_allocate_int32_r2_bounds_type, &
        & pluto_allocate_int32_r2_bounds_name, &
        & pluto_allocate_int32_r2_shape_type, &
        & pluto_allocate_int32_r2_shape_name, &
        & pluto_allocate_label_int32_r2_bounds_type, &
        & pluto_allocate_label_int32_r2_bounds_name, &
        & pluto_allocate_label_int32_r2_shape_type, &
        & pluto_allocate_label_int32_r2_shape_name
    procedure, nopass, private :: pluto_deallocate_int32_r2_type
    procedure, nopass, private :: pluto_deallocate_int32_r2_name
    procedure, nopass, private :: pluto_deallocate_label_int32_r2_type
    procedure, nopass, private :: pluto_deallocate_label_int32_r2_name
    generic :: deallocate => &
        & pluto_deallocate_int32_r2_type, &
        & pluto_deallocate_int32_r2_name, &
        & pluto_deallocate_label_int32_r2_type, &
        & pluto_deallocate_label_int32_r2_name
    procedure, nopass, private :: pluto_allocate_int64_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_int64_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_int64_r2_shape_type
    procedure, nopass, private :: pluto_allocate_int64_r2_shape_name
    procedure, nopass, private :: pluto_allocate_label_int64_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int64_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int64_r2_shape_type
    procedure, nopass, private :: pluto_allocate_label_int64_r2_shape_name
    generic :: allocate => &
        & pluto_allocate_int64_r2_bounds_type, &
        & pluto_allocate_int64_r2_bounds_name, &
        & pluto_allocate_int64_r2_shape_type, &
        & pluto_allocate_int64_r2_shape_name, &
        & pluto_allocate_label_int64_r2_bounds_type, &
        & pluto_allocate_label_int64_r2_bounds_name, &
        & pluto_allocate_label_int64_r2_shape_type, &
        & pluto_allocate_label_int64_r2_shape_name
    procedure, nopass, private :: pluto_deallocate_int64_r2_type
    procedure, nopass, private :: pluto_deallocate_int64_r2_name
    procedure, nopass, private :: pluto_deallocate_label_int64_r2_type
    procedure, nopass, private :: pluto_deallocate_label_int64_r2_name
    generic :: deallocate => &
        & pluto_deallocate_int64_r2_type, &
        & pluto_deallocate_int64_r2_name, &
        & pluto_deallocate_label_int64_r2_type, &
        & pluto_deallocate_label_int64_r2_name
    procedure, nopass, private :: pluto_allocate_real32_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_real32_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_real32_r2_shape_type
    procedure, nopass, private :: pluto_allocate_real32_r2_shape_name
    procedure, nopass, private :: pluto_allocate_label_real32_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real32_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real32_r2_shape_type
    procedure, nopass, private :: pluto_allocate_label_real32_r2_shape_name
    generic :: allocate => &
        & pluto_allocate_real32_r2_bounds_type, &
        & pluto_allocate_real32_r2_bounds_name, &
        & pluto_allocate_real32_r2_shape_type, &
        & pluto_allocate_real32_r2_shape_name, &
        & pluto_allocate_label_real32_r2_bounds_type, &
        & pluto_allocate_label_real32_r2_bounds_name, &
        & pluto_allocate_label_real32_r2_shape_type, &
        & pluto_allocate_label_real32_r2_shape_name
    procedure, nopass, private :: pluto_deallocate_real32_r2_type
    procedure, nopass, private :: pluto_deallocate_real32_r2_name
    procedure, nopass, private :: pluto_deallocate_label_real32_r2_type
    procedure, nopass, private :: pluto_deallocate_label_real32_r2_name
    generic :: deallocate => &
        & pluto_deallocate_real32_r2_type, &
        & pluto_deallocate_real32_r2_name, &
        & pluto_deallocate_label_real32_r2_type, &
        & pluto_deallocate_label_real32_r2_name
    procedure, nopass, private :: pluto_allocate_real64_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_real64_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_real64_r2_shape_type
    procedure, nopass, private :: pluto_allocate_real64_r2_shape_name
    procedure, nopass, private :: pluto_allocate_label_real64_r2_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real64_r2_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real64_r2_shape_type
    procedure, nopass, private :: pluto_allocate_label_real64_r2_shape_name
    generic :: allocate => &
        & pluto_allocate_real64_r2_bounds_type, &
        & pluto_allocate_real64_r2_bounds_name, &
        & pluto_allocate_real64_r2_shape_type, &
        & pluto_allocate_real64_r2_shape_name, &
        & pluto_allocate_label_real64_r2_bounds_type, &
        & pluto_allocate_label_real64_r2_bounds_name, &
        & pluto_allocate_label_real64_r2_shape_type, &
        & pluto_allocate_label_real64_r2_shape_name
    procedure, nopass, private :: pluto_deallocate_real64_r2_type
    procedure, nopass, private :: pluto_deallocate_real64_r2_name
    procedure, nopass, private :: pluto_deallocate_label_real64_r2_type
    procedure, nopass, private :: pluto_deallocate_label_real64_r2_name
    generic :: deallocate => &
        & pluto_deallocate_real64_r2_type, &
        & pluto_deallocate_real64_r2_name, &
        & pluto_deallocate_label_real64_r2_type, &
        & pluto_deallocate_label_real64_r2_name
    procedure, nopass, private :: pluto_allocate_int32_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_int32_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_int32_r3_shape_type
    procedure, nopass, private :: pluto_allocate_int32_r3_shape_name
    procedure, nopass, private :: pluto_allocate_label_int32_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int32_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int32_r3_shape_type
    procedure, nopass, private :: pluto_allocate_label_int32_r3_shape_name
    generic :: allocate => &
        & pluto_allocate_int32_r3_bounds_type, &
        & pluto_allocate_int32_r3_bounds_name, &
        & pluto_allocate_int32_r3_shape_type, &
        & pluto_allocate_int32_r3_shape_name, &
        & pluto_allocate_label_int32_r3_bounds_type, &
        & pluto_allocate_label_int32_r3_bounds_name, &
        & pluto_allocate_label_int32_r3_shape_type, &
        & pluto_allocate_label_int32_r3_shape_name
    procedure, nopass, private :: pluto_deallocate_int32_r3_type
    procedure, nopass, private :: pluto_deallocate_int32_r3_name
    procedure, nopass, private :: pluto_deallocate_label_int32_r3_type
    procedure, nopass, private :: pluto_deallocate_label_int32_r3_name
    generic :: deallocate => &
        & pluto_deallocate_int32_r3_type, &
        & pluto_deallocate_int32_r3_name, &
        & pluto_deallocate_label_int32_r3_type, &
        & pluto_deallocate_label_int32_r3_name
    procedure, nopass, private :: pluto_allocate_int64_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_int64_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_int64_r3_shape_type
    procedure, nopass, private :: pluto_allocate_int64_r3_shape_name
    procedure, nopass, private :: pluto_allocate_label_int64_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int64_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int64_r3_shape_type
    procedure, nopass, private :: pluto_allocate_label_int64_r3_shape_name
    generic :: allocate => &
        & pluto_allocate_int64_r3_bounds_type, &
        & pluto_allocate_int64_r3_bounds_name, &
        & pluto_allocate_int64_r3_shape_type, &
        & pluto_allocate_int64_r3_shape_name, &
        & pluto_allocate_label_int64_r3_bounds_type, &
        & pluto_allocate_label_int64_r3_bounds_name, &
        & pluto_allocate_label_int64_r3_shape_type, &
        & pluto_allocate_label_int64_r3_shape_name
    procedure, nopass, private :: pluto_deallocate_int64_r3_type
    procedure, nopass, private :: pluto_deallocate_int64_r3_name
    procedure, nopass, private :: pluto_deallocate_label_int64_r3_type
    procedure, nopass, private :: pluto_deallocate_label_int64_r3_name
    generic :: deallocate => &
        & pluto_deallocate_int64_r3_type, &
        & pluto_deallocate_int64_r3_name, &
        & pluto_deallocate_label_int64_r3_type, &
        & pluto_deallocate_label_int64_r3_name
    procedure, nopass, private :: pluto_allocate_real32_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_real32_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_real32_r3_shape_type
    procedure, nopass, private :: pluto_allocate_real32_r3_shape_name
    procedure, nopass, private :: pluto_allocate_label_real32_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real32_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real32_r3_shape_type
    procedure, nopass, private :: pluto_allocate_label_real32_r3_shape_name
    generic :: allocate => &
        & pluto_allocate_real32_r3_bounds_type, &
        & pluto_allocate_real32_r3_bounds_name, &
        & pluto_allocate_real32_r3_shape_type, &
        & pluto_allocate_real32_r3_shape_name, &
        & pluto_allocate_label_real32_r3_bounds_type, &
        & pluto_allocate_label_real32_r3_bounds_name, &
        & pluto_allocate_label_real32_r3_shape_type, &
        & pluto_allocate_label_real32_r3_shape_name
    procedure, nopass, private :: pluto_deallocate_real32_r3_type
    procedure, nopass, private :: pluto_deallocate_real32_r3_name
    procedure, nopass, private :: pluto_deallocate_label_real32_r3_type
    procedure, nopass, private :: pluto_deallocate_label_real32_r3_name
    generic :: deallocate => &
        & pluto_deallocate_real32_r3_type, &
        & pluto_deallocate_real32_r3_name, &
        & pluto_deallocate_label_real32_r3_type, &
        & pluto_deallocate_label_real32_r3_name
    procedure, nopass, private :: pluto_allocate_real64_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_real64_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_real64_r3_shape_type
    procedure, nopass, private :: pluto_allocate_real64_r3_shape_name
    procedure, nopass, private :: pluto_allocate_label_real64_r3_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real64_r3_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real64_r3_shape_type
    procedure, nopass, private :: pluto_allocate_label_real64_r3_shape_name
    generic :: allocate => &
        & pluto_allocate_real64_r3_bounds_type, &
        & pluto_allocate_real64_r3_bounds_name, &
        & pluto_allocate_real64_r3_shape_type, &
        & pluto_allocate_real64_r3_shape_name, &
        & pluto_allocate_label_real64_r3_bounds_type, &
        & pluto_allocate_label_real64_r3_bounds_name, &
        & pluto_allocate_label_real64_r3_shape_type, &
        & pluto_allocate_label_real64_r3_shape_name
    procedure, nopass, private :: pluto_deallocate_real64_r3_type
    procedure, nopass, private :: pluto_deallocate_real64_r3_name
    procedure, nopass, private :: pluto_deallocate_label_real64_r3_type
    procedure, nopass, private :: pluto_deallocate_label_real64_r3_name
    generic :: deallocate => &
        & pluto_deallocate_real64_r3_type, &
        & pluto_deallocate_real64_r3_name, &
        & pluto_deallocate_label_real64_r3_type, &
        & pluto_deallocate_label_real64_r3_name
    procedure, nopass, private :: pluto_allocate_int32_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_int32_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_int32_r4_shape_type
    procedure, nopass, private :: pluto_allocate_int32_r4_shape_name
    procedure, nopass, private :: pluto_allocate_label_int32_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int32_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int32_r4_shape_type
    procedure, nopass, private :: pluto_allocate_label_int32_r4_shape_name
    generic :: allocate => &
        & pluto_allocate_int32_r4_bounds_type, &
        & pluto_allocate_int32_r4_bounds_name, &
        & pluto_allocate_int32_r4_shape_type, &
        & pluto_allocate_int32_r4_shape_name, &
        & pluto_allocate_label_int32_r4_bounds_type, &
        & pluto_allocate_label_int32_r4_bounds_name, &
        & pluto_allocate_label_int32_r4_shape_type, &
        & pluto_allocate_label_int32_r4_shape_name
    procedure, nopass, private :: pluto_deallocate_int32_r4_type
    procedure, nopass, private :: pluto_deallocate_int32_r4_name
    procedure, nopass, private :: pluto_deallocate_label_int32_r4_type
    procedure, nopass, private :: pluto_deallocate_label_int32_r4_name
    generic :: deallocate => &
        & pluto_deallocate_int32_r4_type, &
        & pluto_deallocate_int32_r4_name, &
        & pluto_deallocate_label_int32_r4_type, &
        & pluto_deallocate_label_int32_r4_name
    procedure, nopass, private :: pluto_allocate_int64_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_int64_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_int64_r4_shape_type
    procedure, nopass, private :: pluto_allocate_int64_r4_shape_name
    procedure, nopass, private :: pluto_allocate_label_int64_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int64_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int64_r4_shape_type
    procedure, nopass, private :: pluto_allocate_label_int64_r4_shape_name
    generic :: allocate => &
        & pluto_allocate_int64_r4_bounds_type, &
        & pluto_allocate_int64_r4_bounds_name, &
        & pluto_allocate_int64_r4_shape_type, &
        & pluto_allocate_int64_r4_shape_name, &
        & pluto_allocate_label_int64_r4_bounds_type, &
        & pluto_allocate_label_int64_r4_bounds_name, &
        & pluto_allocate_label_int64_r4_shape_type, &
        & pluto_allocate_label_int64_r4_shape_name
    procedure, nopass, private :: pluto_deallocate_int64_r4_type
    procedure, nopass, private :: pluto_deallocate_int64_r4_name
    procedure, nopass, private :: pluto_deallocate_label_int64_r4_type
    procedure, nopass, private :: pluto_deallocate_label_int64_r4_name
    generic :: deallocate => &
        & pluto_deallocate_int64_r4_type, &
        & pluto_deallocate_int64_r4_name, &
        & pluto_deallocate_label_int64_r4_type, &
        & pluto_deallocate_label_int64_r4_name
    procedure, nopass, private :: pluto_allocate_real32_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_real32_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_real32_r4_shape_type
    procedure, nopass, private :: pluto_allocate_real32_r4_shape_name
    procedure, nopass, private :: pluto_allocate_label_real32_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real32_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real32_r4_shape_type
    procedure, nopass, private :: pluto_allocate_label_real32_r4_shape_name
    generic :: allocate => &
        & pluto_allocate_real32_r4_bounds_type, &
        & pluto_allocate_real32_r4_bounds_name, &
        & pluto_allocate_real32_r4_shape_type, &
        & pluto_allocate_real32_r4_shape_name, &
        & pluto_allocate_label_real32_r4_bounds_type, &
        & pluto_allocate_label_real32_r4_bounds_name, &
        & pluto_allocate_label_real32_r4_shape_type, &
        & pluto_allocate_label_real32_r4_shape_name
    procedure, nopass, private :: pluto_deallocate_real32_r4_type
    procedure, nopass, private :: pluto_deallocate_real32_r4_name
    procedure, nopass, private :: pluto_deallocate_label_real32_r4_type
    procedure, nopass, private :: pluto_deallocate_label_real32_r4_name
    generic :: deallocate => &
        & pluto_deallocate_real32_r4_type, &
        & pluto_deallocate_real32_r4_name, &
        & pluto_deallocate_label_real32_r4_type, &
        & pluto_deallocate_label_real32_r4_name
    procedure, nopass, private :: pluto_allocate_real64_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_real64_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_real64_r4_shape_type
    procedure, nopass, private :: pluto_allocate_real64_r4_shape_name
    procedure, nopass, private :: pluto_allocate_label_real64_r4_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real64_r4_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real64_r4_shape_type
    procedure, nopass, private :: pluto_allocate_label_real64_r4_shape_name
    generic :: allocate => &
        & pluto_allocate_real64_r4_bounds_type, &
        & pluto_allocate_real64_r4_bounds_name, &
        & pluto_allocate_real64_r4_shape_type, &
        & pluto_allocate_real64_r4_shape_name, &
        & pluto_allocate_label_real64_r4_bounds_type, &
        & pluto_allocate_label_real64_r4_bounds_name, &
        & pluto_allocate_label_real64_r4_shape_type, &
        & pluto_allocate_label_real64_r4_shape_name
    procedure, nopass, private :: pluto_deallocate_real64_r4_type
    procedure, nopass, private :: pluto_deallocate_real64_r4_name
    procedure, nopass, private :: pluto_deallocate_label_real64_r4_type
    procedure, nopass, private :: pluto_deallocate_label_real64_r4_name
    generic :: deallocate => &
        & pluto_deallocate_real64_r4_type, &
        & pluto_deallocate_real64_r4_name, &
        & pluto_deallocate_label_real64_r4_type, &
        & pluto_deallocate_label_real64_r4_name
    procedure, nopass, private :: pluto_allocate_int32_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_int32_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_int32_r5_shape_type
    procedure, nopass, private :: pluto_allocate_int32_r5_shape_name
    procedure, nopass, private :: pluto_allocate_label_int32_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int32_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int32_r5_shape_type
    procedure, nopass, private :: pluto_allocate_label_int32_r5_shape_name
    generic :: allocate => &
        & pluto_allocate_int32_r5_bounds_type, &
        & pluto_allocate_int32_r5_bounds_name, &
        & pluto_allocate_int32_r5_shape_type, &
        & pluto_allocate_int32_r5_shape_name, &
        & pluto_allocate_label_int32_r5_bounds_type, &
        & pluto_allocate_label_int32_r5_bounds_name, &
        & pluto_allocate_label_int32_r5_shape_type, &
        & pluto_allocate_label_int32_r5_shape_name
    procedure, nopass, private :: pluto_deallocate_int32_r5_type
    procedure, nopass, private :: pluto_deallocate_int32_r5_name
    procedure, nopass, private :: pluto_deallocate_label_int32_r5_type
    procedure, nopass, private :: pluto_deallocate_label_int32_r5_name
    generic :: deallocate => &
        & pluto_deallocate_int32_r5_type, &
        & pluto_deallocate_int32_r5_name, &
        & pluto_deallocate_label_int32_r5_type, &
        & pluto_deallocate_label_int32_r5_name
    procedure, nopass, private :: pluto_allocate_int64_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_int64_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_int64_r5_shape_type
    procedure, nopass, private :: pluto_allocate_int64_r5_shape_name
    procedure, nopass, private :: pluto_allocate_label_int64_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_label_int64_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_label_int64_r5_shape_type
    procedure, nopass, private :: pluto_allocate_label_int64_r5_shape_name
    generic :: allocate => &
        & pluto_allocate_int64_r5_bounds_type, &
        & pluto_allocate_int64_r5_bounds_name, &
        & pluto_allocate_int64_r5_shape_type, &
        & pluto_allocate_int64_r5_shape_name, &
        & pluto_allocate_label_int64_r5_bounds_type, &
        & pluto_allocate_label_int64_r5_bounds_name, &
        & pluto_allocate_label_int64_r5_shape_type, &
        & pluto_allocate_label_int64_r5_shape_name
    procedure, nopass, private :: pluto_deallocate_int64_r5_type
    procedure, nopass, private :: pluto_deallocate_int64_r5_name
    procedure, nopass, private :: pluto_deallocate_label_int64_r5_type
    procedure, nopass, private :: pluto_deallocate_label_int64_r5_name
    generic :: deallocate => &
        & pluto_deallocate_int64_r5_type, &
        & pluto_deallocate_int64_r5_name, &
        & pluto_deallocate_label_int64_r5_type, &
        & pluto_deallocate_label_int64_r5_name
    procedure, nopass, private :: pluto_allocate_real32_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_real32_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_real32_r5_shape_type
    procedure, nopass, private :: pluto_allocate_real32_r5_shape_name
    procedure, nopass, private :: pluto_allocate_label_real32_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real32_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real32_r5_shape_type
    procedure, nopass, private :: pluto_allocate_label_real32_r5_shape_name
    generic :: allocate => &
        & pluto_allocate_real32_r5_bounds_type, &
        & pluto_allocate_real32_r5_bounds_name, &
        & pluto_allocate_real32_r5_shape_type, &
        & pluto_allocate_real32_r5_shape_name, &
        & pluto_allocate_label_real32_r5_bounds_type, &
        & pluto_allocate_label_real32_r5_bounds_name, &
        & pluto_allocate_label_real32_r5_shape_type, &
        & pluto_allocate_label_real32_r5_shape_name
    procedure, nopass, private :: pluto_deallocate_real32_r5_type
    procedure, nopass, private :: pluto_deallocate_real32_r5_name
    procedure, nopass, private :: pluto_deallocate_label_real32_r5_type
    procedure, nopass, private :: pluto_deallocate_label_real32_r5_name
    generic :: deallocate => &
        & pluto_deallocate_real32_r5_type, &
        & pluto_deallocate_real32_r5_name, &
        & pluto_deallocate_label_real32_r5_type, &
        & pluto_deallocate_label_real32_r5_name
    procedure, nopass, private :: pluto_allocate_real64_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_real64_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_real64_r5_shape_type
    procedure, nopass, private :: pluto_allocate_real64_r5_shape_name
    procedure, nopass, private :: pluto_allocate_label_real64_r5_bounds_type
    procedure, nopass, private :: pluto_allocate_label_real64_r5_bounds_name
    procedure, nopass, private :: pluto_allocate_label_real64_r5_shape_type
    procedure, nopass, private :: pluto_allocate_label_real64_r5_shape_name
    generic :: allocate => &
        & pluto_allocate_real64_r5_bounds_type, &
        & pluto_allocate_real64_r5_bounds_name, &
        & pluto_allocate_real64_r5_shape_type, &
        & pluto_allocate_real64_r5_shape_name, &
        & pluto_allocate_label_real64_r5_bounds_type, &
        & pluto_allocate_label_real64_r5_bounds_name, &
        & pluto_allocate_label_real64_r5_shape_type, &
        & pluto_allocate_label_real64_r5_shape_name
    procedure, nopass, private :: pluto_deallocate_real64_r5_type
    procedure, nopass, private :: pluto_deallocate_real64_r5_name
    procedure, nopass, private :: pluto_deallocate_label_real64_r5_type
    procedure, nopass, private :: pluto_deallocate_label_real64_r5_name
    generic :: deallocate => &
        & pluto_deallocate_real64_r5_type, &
        & pluto_deallocate_real64_r5_name, &
        & pluto_deallocate_label_real64_r5_type, &
        & pluto_deallocate_label_real64_r5_name

end type

type(pluto_t) :: pluto

contains

function pluto_devices()
    integer(c_int) :: pluto_devices
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

subroutine pluto_allocate_int32_r1_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int32_r1_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int32_r1_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int32_r1_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r1_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int32_r1_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r1_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r1_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int32_r1_type(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int32_r1_name(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int32_r1_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int32_r1_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int64_r1_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int64_r1_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int64_r1_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int64_r1_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r1_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int64_r1_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r1_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r1_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int64_r1_type(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int64_r1_name(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int64_r1_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int64_r1_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real32_r1_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real32_r1_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real32_r1_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real32_r1_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r1_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real32_r1_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r1_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r1_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real32_r1_type(array, resource)
    real(c_float), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real32_r1_name(array, resource)
    real(c_float), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real32_r1_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real32_r1_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real64_r1_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real64_r1_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real64_r1_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real64_r1_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r1_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real64_r1_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r1_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: shape(1)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r1_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    integer(c_int), intent(in) :: lbounds(1)
    integer(c_int), intent(in) :: ubounds(1)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real64_r1_type(array, resource)
    real(c_double), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real64_r1_name(array, resource)
    real(c_double), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real64_r1_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real64_r1_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int32_r2_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int32_r2_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int32_r2_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int32_r2_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r2_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int32_r2_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r2_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r2_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int32_r2_type(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int32_r2_name(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int32_r2_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int32_r2_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int64_r2_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int64_r2_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int64_r2_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int64_r2_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r2_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int64_r2_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r2_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r2_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int64_r2_type(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int64_r2_name(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int64_r2_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int64_r2_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real32_r2_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real32_r2_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real32_r2_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real32_r2_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r2_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real32_r2_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r2_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r2_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real32_r2_type(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real32_r2_name(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real32_r2_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real32_r2_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real64_r2_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real64_r2_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real64_r2_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real64_r2_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r2_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real64_r2_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r2_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: shape(2)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r2_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    integer(c_int), intent(in) :: lbounds(2)
    integer(c_int), intent(in) :: ubounds(2)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real64_r2_type(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real64_r2_name(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real64_r2_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real64_r2_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int32_r3_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int32_r3_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int32_r3_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int32_r3_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r3_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int32_r3_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r3_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r3_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int32_r3_type(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int32_r3_name(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int32_r3_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int32_r3_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int64_r3_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int64_r3_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int64_r3_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int64_r3_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r3_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int64_r3_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r3_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r3_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int64_r3_type(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int64_r3_name(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int64_r3_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int64_r3_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real32_r3_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real32_r3_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real32_r3_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real32_r3_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r3_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real32_r3_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r3_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r3_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real32_r3_type(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real32_r3_name(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real32_r3_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real32_r3_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real64_r3_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real64_r3_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real64_r3_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real64_r3_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r3_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real64_r3_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r3_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: shape(3)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r3_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    integer(c_int), intent(in) :: lbounds(3)
    integer(c_int), intent(in) :: ubounds(3)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real64_r3_type(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real64_r3_name(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real64_r3_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real64_r3_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int32_r4_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int32_r4_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int32_r4_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int32_r4_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r4_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int32_r4_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r4_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r4_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int32_r4_type(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int32_r4_name(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int32_r4_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int32_r4_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int64_r4_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int64_r4_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int64_r4_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int64_r4_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r4_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int64_r4_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r4_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r4_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int64_r4_type(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int64_r4_name(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int64_r4_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int64_r4_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real32_r4_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real32_r4_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real32_r4_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real32_r4_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r4_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real32_r4_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r4_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r4_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real32_r4_type(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real32_r4_name(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real32_r4_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real32_r4_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real64_r4_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real64_r4_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real64_r4_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real64_r4_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r4_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real64_r4_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r4_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: shape(4)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r4_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    integer(c_int), intent(in) :: lbounds(4)
    integer(c_int), intent(in) :: ubounds(4)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real64_r4_type(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real64_r4_name(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real64_r4_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real64_r4_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int32_r5_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int32_r5_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int32_r5_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int32_r5_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r5_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int32_r5_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int32_r5_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int32_r5_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int32_r5_type(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int32_r5_name(array, resource)
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int32_r5_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int32_r5_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int32_t), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_int64_r5_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_int64_r5_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_int64_r5_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_int64_r5_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r5_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_int64_r5_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_int64_r5_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_int64_r5_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_int64_r5_type(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_int64_r5_name(array, resource)
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_int64_r5_type(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_int64_r5_name(label, array, resource)
    character(len=*), intent(in) :: label
    integer(c_int64_t), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real32_r5_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real32_r5_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real32_r5_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real32_r5_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r5_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real32_r5_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real32_r5_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real32_r5_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real32_r5_type(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real32_r5_name(array, resource)
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real32_r5_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real32_r5_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_float), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine
subroutine pluto_allocate_real64_r5_bounds_type(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_real64_r5_bounds_name(array, lbounds, ubounds, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_real64_r5_shape_type(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(array, shape, resource)
end subroutine

subroutine pluto_allocate_real64_r5_shape_name(array, shape, resource)
    use pluto_module_allocate_deallocate, only : pluto_allocate
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r5_bounds_type(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, resource)
end subroutine

subroutine pluto_allocate_label_real64_r5_shape_name(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, shape, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_allocate_label_real64_r5_shape_type(label, array, shape, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: shape(5)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_allocate(label, array, shape, resource)
end subroutine

subroutine pluto_allocate_label_real64_r5_bounds_name(label, array, lbounds, ubounds, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    integer(c_int), intent(in) :: lbounds(5)
    integer(c_int), intent(in) :: ubounds(5)
    character(len=*), intent(in) :: resource
    call pluto_allocate(label, array, lbounds, ubounds, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_real64_r5_type(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(array, resource)
end subroutine

subroutine pluto_deallocate_real64_r5_name(array, resource)
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(array, pluto_get_registered_resource(resource))
end subroutine

subroutine pluto_deallocate_label_real64_r5_type(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    type(pluto_memory_resource), intent(in) :: resource
    call pluto_deallocate(label, array, resource)
end subroutine

subroutine pluto_deallocate_label_real64_r5_name(label, array, resource)
    character(len=*), intent(in) :: label
    real(c_double), pointer, intent(inout) :: array(:,:,:,:,:)
    character(len=*), intent(in) :: resource
    call pluto_deallocate(label, array, pluto_get_registered_resource(resource))
end subroutine

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

end module
