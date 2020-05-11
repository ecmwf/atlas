! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_metadata_module

use fckit_object_module, only : fckit_object

implicit none

private :: fckit_object

public :: atlas_Metadata

private

integer, parameter :: MAX_STR_LEN=255

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_Metadata

! Purpose :
! -------
!   *Metadata* : Container of Metadata, parameters or attributes
!       The Metadata are seted as key, value pairs

! Methods :
! -------
!   set : set a new property with given key and value
!   set : Modify a property with given key and value
!   get : Return a property value for given key

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: set_logical => Metadata__set_logical
  procedure, private :: set_int32 => Metadata__set_int32
  procedure, private :: set_real32 => Metadata__set_real32
  procedure, private :: set_real64 => Metadata__set_real64
  procedure, private :: set_string => Metadata__set_string
  procedure, private :: set_array_int32 => Metadata__set_array_int32
  procedure, private :: set_array_int64 => Metadata__set_array_int64
  procedure, private :: set_array_real32 => Metadata__set_array_real32
  procedure, private :: set_array_real64 => Metadata__set_array_real64
  procedure :: has => Metadata__has
  generic :: set => set_logical, set_int32, set_real32, set_real64, set_string, set_array_int32, &
    set_array_int64, set_array_real32, set_array_real64
  procedure :: get_int32 => Metadata__get_int32
  procedure :: get_logical => Metadata__get_logical
  procedure :: get_real32 => Metadata__get_real32
  procedure :: get_real64 => Metadata__get_real64
  procedure :: get_string => Metadata__get_string
  procedure :: get_array_int32 => Metadata__get_array_int32
  procedure :: get_array_int64 => Metadata__get_array_int64
  procedure :: get_array_real32 => Metadata__get_array_real32
  procedure :: get_array_real64 => Metadata__get_array_real64
  generic :: get => get_int32, get_logical, get_real32, get_real64, get_string, get_array_int32, &
    get_array_int64, get_array_real32, get_array_real64

  procedure :: print => Metadata__print
  procedure :: json => Metadata__json

  procedure, public :: delete => atlas_Metadata__delete

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Metadata__final_auto
#endif

END TYPE atlas_Metadata

!------------------------------------------------------------------------------

interface atlas_Metadata
  module procedure atlas_Metadata__ctor
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! Metadata routines

function atlas_Metadata__ctor() result(this)
  use atlas_metadata_c_binding
  use fckit_c_interop_module
  type(atlas_Metadata) :: this
  call this%reset_c_ptr( atlas__Metadata__new() )
end function atlas_Metadata__ctor

subroutine atlas_Metadata__delete(this)
  use atlas_metadata_c_binding
  class(atlas_Metadata), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Metadata__delete(this%CPTR_PGIBUG_A)
  end if
  call this%reset_c_ptr()
end subroutine atlas_Metadata__delete

function Metadata__has(this, name) result(value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical :: value
  integer :: value_int
  value_int =  atlas__Metadata__has(this%CPTR_PGIBUG_A, c_str(name) )
  if( value_int == 1 ) then
    value = .True.
  else
    value = .False.
  end if
end function Metadata__has

subroutine Metadata__set_logical(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical, intent(in) :: value
  integer :: value_int
  if( value ) then
    value_int = 1
  else
    value_int = 0
  end if
  call atlas__Metadata__set_int(this%CPTR_PGIBUG_A, c_str(name), value_int )
end subroutine Metadata__set_logical

subroutine Metadata__set_int32(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call atlas__Metadata__set_int(this%CPTR_PGIBUG_A, c_str(name), value)
end subroutine Metadata__set_int32

subroutine Metadata__set_real32(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Metadata__set_float(this%CPTR_PGIBUG_A, c_str(name) ,value)
end subroutine Metadata__set_real32

subroutine Metadata__set_real64(this, name, value)
  use atlas_metadata_c_binding
  use, intrinsic :: iso_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Metadata__set_double(this%CPTR_PGIBUG_A, c_str(name) ,value)
end subroutine Metadata__set_real64

subroutine Metadata__set_string(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Metadata__set_string(this%CPTR_PGIBUG_A, c_str(name) , c_str(value) )
end subroutine Metadata__set_string

subroutine Metadata__get_logical(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(out) :: value
  integer :: value_int
  value_int = atlas__Metadata__get_int(this%CPTR_PGIBUG_A,c_str(name) )
  if (value_int > 0) then
    value = .True.
  else
    value = .False.
  end if
end subroutine Metadata__get_logical

subroutine Metadata__get_int32(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(out) :: value
  value = atlas__Metadata__get_int(this%CPTR_PGIBUG_A, c_str(name) )
end subroutine Metadata__get_int32

subroutine Metadata__get_real32(this, name, value)
  use atlas_metadata_c_binding
  use, intrinsic :: iso_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(out) :: value
  value = atlas__Metadata__get_float(this%CPTR_PGIBUG_A, c_str(name) )
end subroutine Metadata__get_real32

subroutine Metadata__get_real64(this, name, value)
  use atlas_metadata_c_binding
  use, intrinsic :: iso_c_binding
  use fckit_c_interop_module, only : c_str
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(out) :: value
  value = atlas__Metadata__get_double(this%CPTR_PGIBUG_A, c_str(name) )
end subroutine Metadata__get_real64

subroutine Metadata__get_string(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str, c_str_to_string
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(out) :: value
  character(len=MAX_STR_LEN) :: value_cstr
  call atlas__Metadata__get_string(this%CPTR_PGIBUG_A, c_str(name), value_cstr, MAX_STR_LEN )
  value = c_str_to_string(value_cstr)
end subroutine Metadata__get_string

subroutine Metadata__set_array_int32(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: value(:)
  call atlas__Metadata__set_array_int(this%CPTR_PGIBUG_A, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_int32

subroutine Metadata__set_array_int64(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: value(:)
  call atlas__Metadata__set_array_long(this%CPTR_PGIBUG_A, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_int64

subroutine Metadata__set_array_real32(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value(:)
  call atlas__Metadata__set_array_float(this%CPTR_PGIBUG_A, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_real32

subroutine Metadata__set_array_real64(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value(:)
  call atlas__Metadata__set_array_double(this%CPTR_PGIBUG_A, c_str(name), &
    & value, size(value) )
end subroutine Metadata__set_array_real64

subroutine Metadata__get_array_int32(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_int), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_int(this%CPTR_PGIBUG_A, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,[value_size])
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call c_ptr_free(value_cptr)
end subroutine Metadata__get_array_int32

subroutine Metadata__get_array_int64(this, name, value)
  use atlas_metadata_c_binding
  use fckit_c_interop_module, only : c_str, c_ptr_free
  use, intrinsic :: iso_c_binding
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  integer(c_long), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_long(this%CPTR_PGIBUG_A, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call c_ptr_free(value_cptr)
end subroutine Metadata__get_array_int64

subroutine Metadata__get_array_real32(this, name, value)
  use atlas_metadata_c_binding
  use, intrinsic :: iso_c_binding
  use fckit_c_interop_module, only : c_str, c_ptr_free
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_float), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_float(this%CPTR_PGIBUG_A, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call c_ptr_free(value_cptr)
end subroutine Metadata__get_array_real32

subroutine Metadata__get_array_real64(this, name, value)
  use atlas_metadata_c_binding
  use, intrinsic :: iso_c_binding
  use fckit_c_interop_module, only : c_str, c_ptr_free
  class(atlas_Metadata), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), allocatable, intent(out) :: value(:)
  type(c_ptr) :: value_cptr
  real(c_double), pointer :: value_fptr(:)
  integer :: value_size
  integer :: value_allocated
  call atlas__Metadata__get_array_double(this%CPTR_PGIBUG_A, c_str(name), &
    & value_cptr, value_size, value_allocated )
  call c_f_pointer(value_cptr,value_fptr,(/value_size/))
  allocate(value(value_size))
  value(:) = value_fptr(:)
  if( value_allocated == 1 ) call c_ptr_free(value_cptr)
end subroutine Metadata__get_array_real64

subroutine MetaData__print(this,channel)
  use atlas_metadata_c_binding
  use fckit_log_module, only : fckit_logchannel
  class(atlas_Metadata), intent(in) :: this
  type(fckit_logchannel), intent(in) :: channel
  call atlas__Metadata__print(this%CPTR_PGIBUG_A,channel%CPTR_PGIBUG_A)
end subroutine Metadata__print

function Metadata__json(this) result(json)
  use atlas_metadata_c_binding
  use, intrinsic :: iso_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  character(len=:), allocatable :: json
  class(atlas_Metadata), intent(in) :: this
  type(c_ptr) :: json_cptr
  integer(c_int) :: json_size
  integer(c_int) :: json_allocated
  call atlas__Metadata__json(this%CPTR_PGIBUG_A,json_cptr,json_size,json_allocated)
  allocate(character(len=json_size) :: json )
  json = c_ptr_to_string(json_cptr)
  if( json_allocated == 1 ) call c_ptr_free(json_cptr)
end function Metadata__json

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Metadata__final_auto(this)
  type(atlas_Metadata), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Metadata__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_metadata_module

