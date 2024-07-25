! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_config_module

use fckit_shared_object_module, only : fckit_shared_object, fckit_c_deleter, fckit_c_nodeleter

implicit none

public :: atlas_Config

private

TYPE, extends(fckit_shared_object) :: atlas_Config

! Purpose :
! -------
!   *Config* : Container of Config, parameters or attributes
!       The Config are seted as key, value pairs

! Methods :
! -------
!   set : set a new property with given key and value
!   set : Modify a property with given key and value
!   get : Return a property value for given key

! Author :
! ------
!   June-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: set_config => atlas_Config__set_config
  procedure, private :: set_config_list => atlas_Config__set_config_list
  procedure, private :: set_copy => atlas_Config__copy
  procedure, private :: set_logical => atlas_Config__set_logical
  procedure, private :: set_int32 => atlas_Config__set_int32
  procedure, private :: set_int64 => atlas_Config__set_int64
  procedure, private :: set_real32 => atlas_Config__set_real32
  procedure, private :: set_real64 => atlas_Config__set_real64
  procedure, private :: set_string => atlas_Config__set_string
  procedure, private :: set_array_int32 => atlas_Config__set_array_int32
  procedure, private :: set_array_int64 => atlas_Config__set_array_int64
  procedure, private :: set_array_real32 => atlas_Config__set_array_real32
  procedure, private :: set_array_real64 => atlas_Config__set_array_real64
  procedure :: has => atlas_Config__has
  generic :: set => set_config, set_config_list, set_logical, set_int32, set_int64, set_real32, set_real64, &
                    set_string, set_array_int32, set_array_int64, set_array_real32, set_array_real64, set_copy
  procedure, private :: get_config => atlas_Config__get_config
  procedure, private :: get_config_list => atlas_Config__get_config_list
  procedure, private :: get_int32 => atlas_Config__get_int32
  procedure, private :: get_logical => atlas_Config__get_logical
  procedure, private :: get_real32 => atlas_Config__get_real32
  procedure, private :: get_real64 => atlas_Config__get_real64
  procedure, private :: get_string => atlas_Config__get_string
  procedure, private :: get_array_int32 => atlas_Config__get_array_int32
  procedure, private :: get_array_int64 => atlas_Config__get_array_int64
  procedure, private :: get_array_real32 => atlas_Config__get_array_real32
  procedure, private :: get_array_real64 => atlas_Config__get_array_real64
  generic :: get => get_config, get_config_list, get_int32, get_logical, get_real32, get_real64, &
                    get_string, get_array_int32, get_array_int64, get_array_real32, get_array_real64
  procedure :: json => atlas_Config__json

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Config__final_auto
#endif

END TYPE atlas_Config

!------------------------------------------------------------------------------

interface atlas_Config
  module procedure atlas_Config__ctor
  module procedure atlas_Config__ctor_from_cptr
  module procedure atlas_Config__ctor_from_file
  module procedure atlas_Config__ctor_from_json
end interface

!------------------------------------------------------------------------------

private :: fckit_shared_object
private :: fckit_c_deleter
private :: fckit_c_nodeleter

!========================================================
contains
!========================================================

! -----------------------------------------------------------------------------
! Config routines

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Config__final_auto(this)
  type(atlas_Config), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Config__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

function atlas_Config__ctor() result(this)
  use atlas_Config_c_binding
  type(atlas_Config) :: this
  call this%reset_c_ptr( atlas__Config__new(), &
    & fckit_c_deleter(atlas__Config__delete) )
  call this%return()
end function atlas_Config__ctor

function atlas_Config__ctor_from_cptr(cptr,delete) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_Config_c_binding
  type(c_ptr), value :: cptr
  type(atlas_Config) :: this
  logical, optional :: delete
  logical :: opt_delete
  opt_delete = .true.
  if( present(delete) ) opt_delete = delete
  if( opt_delete ) then
    call this%reset_c_ptr( cptr, fckit_c_deleter(atlas__Config__delete) )
  else
    call this%reset_c_ptr( cptr, fckit_c_nodeleter() )
  endif
  call this%return()
end function

function atlas_Config__ctor_from_json(json) result(this)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  use atlas_JSON_module
  type(atlas_Config) :: this
  class(atlas_JSON) :: json
  call this%reset_c_ptr( atlas__Config__new_from_json(c_str(json%str())), &
    & fckit_c_deleter(atlas__Config__delete) )
  call this%return()
end function atlas_Config__ctor_from_json

function atlas_Config__ctor_from_file(path) result(this)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  use atlas_JSON_module
  type(atlas_Config) :: this
  class(atlas_PathName), intent(in) :: path
  call this%reset_c_ptr( atlas__Config__new_from_file(c_str(path%str())), &
    & fckit_c_deleter(atlas__Config__delete) )
  call this%return()
end function atlas_Config__ctor_from_file

function atlas_Config__has(this, name) result(value)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical :: value
  integer :: value_int
  value_int =  atlas__Config__has(this%CPTR_PGIBUG_B, c_str(name) )
  if( value_int == 1 ) then
    value = .True.
  else
    value = .False.
  end if
end function atlas_Config__has

subroutine atlas_Config__copy(this, config)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  class(atlas_Config), intent(in) :: config
  call atlas__Config__copy(this%CPTR_PGIBUG_B, config%CPTR_PGIBUG_B)
end subroutine atlas_Config__copy

subroutine atlas_Config__set_config(this, name, value)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  class(atlas_Config), intent(in) :: value
  call atlas__Config__set_config(this%CPTR_PGIBUG_B, c_str(name), value%CPTR_PGIBUG_B )
end subroutine atlas_Config__set_config

subroutine atlas_Config__set_config_list(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_loc
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(atlas_Config), intent(in) :: value(:) !PGI (17.7) compiler bug when "type" replaced with "class"
  type(c_ptr), target :: value_cptrs(size(value))
  integer :: j
  if( size(value) > 0 ) then
    do j=1,size(value)
      value_cptrs(j) = value(j)%CPTR_PGIBUG_B
    enddo
    call atlas__Config__set_config_list(this%CPTR_PGIBUG_B, c_str(name), c_loc(value_cptrs(1)), size(value_cptrs) )
  endif
end subroutine atlas_Config__set_config_list


subroutine atlas_Config__set_logical(this, name, value)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical, intent(in) :: value
  integer :: value_int
  if( value ) then
    value_int = 1
  else
    value_int = 0
  end if
  call atlas__Config__set_int(this%CPTR_PGIBUG_B, c_str(name), value_int )
end subroutine atlas_Config__set_logical

subroutine atlas_Config__set_int32(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_int
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: value
  call atlas__Config__set_int(this%CPTR_PGIBUG_B, c_str(name), value)
end subroutine atlas_Config__set_int32

subroutine atlas_Config__set_int64(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_long
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: value
  call atlas__Config__set_long(this%CPTR_PGIBUG_B, c_str(name), value)
end subroutine atlas_Config__set_int64

subroutine atlas_Config__set_real32(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_float
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Config__set_float(this%CPTR_PGIBUG_B, c_str(name) ,value)
end subroutine atlas_Config__set_real32

subroutine atlas_Config__set_real64(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Config__set_double(this%CPTR_PGIBUG_B, c_str(name) ,value)
end subroutine atlas_Config__set_real64

subroutine atlas_Config__set_string(this, name, value)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Config__set_string(this%CPTR_PGIBUG_B, c_str(name) , c_str(value) )
end subroutine atlas_Config__set_string

function atlas_Config__get_config(this, name, value) result(found)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  class(atlas_Config), intent(inout) :: value
  integer :: found_int
  value = atlas_Config()
  found_int = atlas__Config__get_config( this%CPTR_PGIBUG_B, c_str(name), value%CPTR_PGIBUG_B )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_config

function atlas_Config__get_config_list(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer, c_null_ptr
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Config), allocatable, intent(inout) :: value(:)
  type(c_ptr) :: value_list_cptr
  type(c_ptr), pointer :: value_cptrs(:)
  integer :: value_list_allocated
  integer :: value_list_size
  integer :: found_int
  integer :: j
  value_list_cptr = c_null_ptr
  found_int = atlas__Config__get_config_list(this%CPTR_PGIBUG_B, c_str(name), &
    & value_list_cptr, value_list_size, value_list_allocated )
  if( found_int == 1 ) then
    call c_f_pointer(value_list_cptr,value_cptrs,(/value_list_size/))
    if( allocated(value) ) deallocate(value)
    allocate(value(value_list_size))
    do j=1,value_list_size
      value(j) = atlas_Config( value_cptrs(j) )
    enddo
    if( value_list_allocated == 1 ) call c_ptr_free(value_list_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_config_list

function atlas_Config__get_logical(this, name, value) result(found)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(inout) :: value
  integer :: value_int
  integer :: found_int
  found_int = atlas__Config__get_int(this%CPTR_PGIBUG_B,c_str(name), value_int )
  found = .False.
  if (found_int == 1) found = .True.
  if (found) then
    if (value_int > 0) then
      value = .True.
    else
      value = .False.
    end if
  endif
end function atlas_Config__get_logical

function atlas_Config__get_int32(this, name, value) result(found)
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(inout) :: value
  integer :: found_int
  found_int = atlas__Config__get_int(this%CPTR_PGIBUG_B, c_str(name), value )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_int32

function atlas_Config__get_real32(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_float
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(inout) :: value
  integer :: found_int
  found_int = atlas__Config__get_float(this%CPTR_PGIBUG_B, c_str(name), value )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_real32

function atlas_Config__get_real64(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(inout) :: value
  integer :: found_int
  found_int = atlas__Config__get_double(this%CPTR_PGIBUG_B, c_str(name), value )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_real64

function atlas_Config__get_string(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(inout) :: value
  type(c_ptr) :: value_cptr
  integer :: found_int
  integer(c_int) :: value_size
  integer(c_int) :: value_allocated
  found_int = atlas__Config__get_string(this%CPTR_PGIBUG_B,c_str(name),value_cptr,value_size,value_allocated)
  if( found_int == 1 ) then
    if( allocated(value) ) deallocate(value)
    allocate(character(len=value_size) :: value )
    value = c_ptr_to_string(value_cptr)
    if( value_allocated == 1 ) call c_ptr_free(value_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_string

subroutine atlas_Config__set_array_int32(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_int
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: value(:)
  call atlas__Config__set_array_int(this%CPTR_PGIBUG_B, c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_int32

subroutine atlas_Config__set_array_int64(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_long
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: value(:)
  call atlas__Config__set_array_long(this%CPTR_PGIBUG_B, c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_int64

subroutine atlas_Config__set_array_real32(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_float
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value(:)
  call atlas__Config__set_array_float(this%CPTR_PGIBUG_B, c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_real32

subroutine atlas_Config__set_array_real64(this, name, value)
  use, intrinsic :: iso_c_binding, only : c_double
  use fckit_c_interop_module, only : c_str
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value(:)
  call atlas__Config__set_array_double(this%CPTR_PGIBUG_B, c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_real64

function atlas_Config__get_array_int32(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_f_pointer
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), allocatable, intent(inout) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_int), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  integer :: found_int
  found_int = atlas__Config__get_array_int(this%CPTR_PGIBUG_B, c_str(name), &
    & value_cptr, value_size, value_allocated )
  if (found_int ==1 ) then
    call c_f_pointer(value_cptr,value_fptr,(/value_size/))
    if( allocated(value) ) deallocate(value)
    allocate(value(value_size))
    value(:) = value_fptr(:)
    if( value_allocated == 1 ) call c_ptr_free(value_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_array_int32

function atlas_Config__get_array_int64(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_long, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), allocatable, intent(inout) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_long), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  integer :: found_int
  found_int = atlas__Config__get_array_long(this%CPTR_PGIBUG_B, c_str(name), &
    & value_cptr, value_size, value_allocated )
  if (found_int == 1) then
    call c_f_pointer(value_cptr,value_fptr,(/value_size/))
    if( allocated(value) ) deallocate(value)
    allocate(value(value_size))
    value(:) = value_fptr(:)
    if( value_allocated == 1 ) call c_ptr_free(value_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_array_int64

function atlas_Config__get_array_real32(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_float, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), allocatable, intent(inout) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_float), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  integer :: found_int
  found_int = atlas__Config__get_array_float(this%CPTR_PGIBUG_B, c_str(name), &
    & value_cptr, value_size, value_allocated )
  if (found_int == 1 ) then
    call c_f_pointer(value_cptr,value_fptr,(/value_size/))
    if( allocated(value) ) deallocate(value)
    allocate(value(value_size))
    value(:) = value_fptr(:)
    if( value_allocated == 1 ) call c_ptr_free(value_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_array_real32

function atlas_Config__get_array_real64(this, name, value) result(found)
  use, intrinsic :: iso_c_binding, only : c_double, c_ptr, c_f_pointer
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), allocatable, intent(inout) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_double), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  integer :: found_int
  found_int = atlas__Config__get_array_double(this%CPTR_PGIBUG_B, c_str(name), &
    & value_cptr, value_size, value_allocated )
  if (found_int == 1) then
    call c_f_pointer(value_cptr,value_fptr,(/value_size/))
    if( allocated(value) ) deallocate(value)
    allocate(value(value_size))
    value(:) = value_fptr(:)
    if( value_allocated == 1 ) call c_ptr_free(value_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_array_real64

function atlas_Config__json(this) result(json)
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
  use atlas_Config_c_binding
  character(len=:), allocatable :: json
  class(atlas_Config), intent(in) :: this
  type(c_ptr) :: json_cptr
  integer(c_int) :: json_size
  integer(c_int) :: json_allocated
  call atlas__Config__json(this%CPTR_PGIBUG_B,json_cptr,json_size,json_allocated)
  allocate(character(len=json_size) :: json )
  json = c_ptr_to_string(json_cptr)
  if( json_allocated == 1 ) call c_ptr_free(json_cptr)
end function atlas_Config__json

end module atlas_config_module

