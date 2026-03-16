! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_host

use, intrinsic :: iso_c_binding, only : c_int, c_double, c_float, c_int32_t, c_int64_t, c_loc
use pluto_module_allocate_deallocate, only : pluto_allocate, pluto_deallocate
use pluto_module_allocator, only : pluto_allocator
use pluto_module_memory_resource, only : pluto_memory_resource

implicit none
private

! A workaround for NVHPC compiler bug.
! If the compiler is NVHPC we need the THIS argument to be passed explicitly.
! There are otherwise issues with ambiguous generic interfaces between pluto_host_t and pluto_device_t,
! which are resolved by passing the THIS argument explicitly.

#if !defined(__NVCOMPILER)
#define NOPASS , nopass
#define THIS
#define THIS_COMMA
#define CLASS_THIS
#else
#define NOPASS
#define THIS this
#define THIS_COMMA this,
#define CLASS_THIS class(pluto_host_t), intent(in) :: this
#endif

public :: pluto_host_t

type pluto_host_t
   integer :: dummy
contains
    procedure NOPASS :: get_default_resource  => pluto_host_get_default_resource
    procedure NOPASS :: make_allocator        => pluto_host_make_allocator
    procedure NOPASS, private :: pluto_host_set_default_resource_name
    procedure NOPASS, private :: pluto_host_set_default_resource_type
    generic, public :: set_default_resource => pluto_host_set_default_resource_type, &
                                             & pluto_host_set_default_resource_name
end type

contains

function pluto_host_get_default_resource(THIS) result(memory_resource)
    type(pluto_memory_resource) :: memory_resource
    CLASS_THIS
    interface
        function c_pluto_host_get_default_resource() result(mr) bind(c)
            use iso_c_binding, only: c_ptr
            type(c_ptr) :: mr 
        end function
    end interface
    memory_resource%c_memory_resource = c_pluto_host_get_default_resource()
end function
subroutine pluto_host_set_default_resource_name(THIS_COMMA name)
    CLASS_THIS
    character(len=*), target, intent(in) :: name
    interface
        subroutine c_pluto_host_set_default_resource_name(name, name_size) bind(c)
            use iso_c_binding, only: c_ptr, c_int
            type(c_ptr), value, intent(in) :: name
            integer(c_int), value, intent(in) :: name_size
        end subroutine
    end interface
    call c_pluto_host_set_default_resource_name(c_loc(name), len_trim(name,kind=c_int))
end subroutine

subroutine pluto_host_set_default_resource_type(THIS_COMMA memory_resource)
    CLASS_THIS
    type(pluto_memory_resource), intent(in) :: memory_resource
    interface
        subroutine c_pluto_host_set_default_resource_ptr(mr) bind(c)
            use iso_c_binding, only: c_ptr
            type(c_ptr), value, intent(in) :: mr
        end subroutine
    end interface
    call c_pluto_host_set_default_resource_ptr(memory_resource%c_memory_resource)
end subroutine


function pluto_host_make_allocator(THIS) result(allocator)
    type(pluto_allocator) :: allocator
    CLASS_THIS
    allocator%memory_resource = pluto_host_get_default_resource(THIS)
end function

end module
