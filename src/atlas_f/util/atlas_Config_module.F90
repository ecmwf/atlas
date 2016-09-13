
module atlas_config_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_long, c_double, c_float, c_f_pointer, c_null_ptr, c_loc
use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use fckit_refcounted_module, only : fckit_refcounted_fortran

implicit none

private :: fckit_refcounted_fortran
private :: c_ptr, c_int, c_long, c_double, c_float, c_f_pointer, c_null_ptr, c_loc
private :: c_str, c_ptr_to_string, c_ptr_free
public :: atlas_Config

private

TYPE, extends(fckit_refcounted_fortran) :: atlas_Config

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
  procedure, private :: set_logical => atlas_Config__set_logical
  procedure, private :: set_int32 => atlas_Config__set_int32
  procedure, private :: set_real32 => atlas_Config__set_real32
  procedure, private :: set_real64 => atlas_Config__set_real64
  procedure, private :: set_string => atlas_Config__set_string
  procedure, private :: set_array_int32 => atlas_Config__set_array_int32
  procedure, private :: set_array_int64 => atlas_Config__set_array_int64
  procedure, private :: set_array_real32 => atlas_Config__set_array_real32
  procedure, private :: set_array_real64 => atlas_Config__set_array_real64
  procedure :: has => atlas_Config__has
  generic :: set => set_config, set_config_list, set_logical, set_int32, set_real32, set_real64, &
                    set_string, set_array_int32, set_array_int64, set_array_real32, set_array_real64
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

  procedure, public :: delete => atlas_Config__delete

END TYPE atlas_Config

!------------------------------------------------------------------------------

interface atlas_Config
  module procedure atlas_Config__ctor
  module procedure atlas_Config__ctor_from_file
  module procedure atlas_Config__ctor_from_json
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! Config routines

function atlas_Config__ctor() result(Config)
  use atlas_Config_c_binding
  type(atlas_Config) :: Config
  call Config%reset_c_ptr( atlas__Config__new() )
end function atlas_Config__ctor

function atlas_Config__ctor_from_json(json) result(Config)
  use atlas_Config_c_binding
  use atlas_JSON_module
  type(atlas_Config) :: Config
  class(atlas_JSON) :: json
  call Config%reset_c_ptr( atlas__Config__new_from_json(c_str(json%str())) )
end function atlas_Config__ctor_from_json

function atlas_Config__ctor_from_file(path) result(Config)
  use atlas_Config_c_binding
  use atlas_JSON_module
  type(atlas_Config) :: Config
  class(atlas_PathName), intent(in) :: path
  call Config%reset_c_ptr( atlas__Config__new_from_file(c_str(path%str())) )
end function atlas_Config__ctor_from_file

subroutine atlas_Config__delete(this)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Config__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_Config__delete


function atlas_Config__has(this, name) result(value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical :: value
  integer :: value_int
  value_int =  atlas__Config__has(this%c_ptr(), c_str(name) )
  if( value_int == 1 ) then
    value = .True.
  else
    value = .False.
  end if
end function atlas_Config__has

subroutine atlas_Config__set_config(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  class(atlas_Config), intent(in) :: value
  call atlas__Config__set_config(this%c_ptr(), c_str(name), value%c_ptr() )
end subroutine atlas_Config__set_config

subroutine atlas_Config__set_config_list(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  class(atlas_Config), intent(in) :: value(:)
  type(c_ptr), target :: value_cptrs(size(value))
  integer :: j
  if( size(value) > 0 ) then
    do j=1,size(value)
      value_cptrs(j) = value(j)%c_ptr()
    enddo
    call atlas__Config__set_config_list(this%c_ptr(), c_str(name), c_loc(value_cptrs(1)), size(value_cptrs) )
  endif
end subroutine atlas_Config__set_config_list

subroutine atlas_Config__set_logical(this, name, value)
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
  call atlas__Config__set_int(this%c_ptr(), c_str(name), value_int )
end subroutine atlas_Config__set_logical

subroutine atlas_Config__set_int32(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call atlas__Config__set_int(this%c_ptr(), c_str(name), value)
end subroutine atlas_Config__set_int32

subroutine atlas_Config__set_real32(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Config__set_float(this%c_ptr(), c_str(name) ,value)
end subroutine atlas_Config__set_real32

subroutine atlas_Config__set_real64(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Config__set_double(this%c_ptr(), c_str(name) ,value)
end subroutine atlas_Config__set_real64

subroutine atlas_Config__set_string(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Config__set_string(this%c_ptr(), c_str(name) , c_str(value) )
end subroutine atlas_Config__set_string

function atlas_Config__get_config(this, name, value) result(found)
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  class(atlas_Config), intent(inout) :: value
  integer :: found_int
  if( value%is_null() ) then
    call value%reset_c_ptr( atlas__Config__new() )
  endif
  found_int = atlas__Config__get_config(this%c_ptr(), c_str(name), value%c_ptr() )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_config

function atlas_Config__get_config_list(this, name, value) result(found)
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
  found_int = atlas__Config__get_config_list(this%c_ptr(), c_str(name), &
    & value_list_cptr, value_list_size, value_list_allocated )
  if( found_int == 1 ) then
    call c_f_pointer(value_list_cptr,value_cptrs,(/value_list_size/))
    if( allocated(value) ) deallocate(value)
    allocate(value(value_list_size))
    do j=1,value_list_size
      call value(j)%reset_c_ptr( value_cptrs(j) )
    enddo
    if( value_list_allocated == 1 ) call c_ptr_free(value_list_cptr)
  endif
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_config_list

function atlas_Config__get_logical(this, name, value) result(found)
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(inout) :: value
  integer :: value_int
  integer :: found_int
  found_int = atlas__Config__get_int(this%c_ptr(),c_str(name), value_int )
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
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(inout) :: value
  integer :: found_int
  found_int = atlas__Config__get_int(this%c_ptr(), c_str(name), value )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_int32

function atlas_Config__get_real32(this, name, value) result(found)
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(inout) :: value
  integer :: found_int
  found_int = atlas__Config__get_float(this%c_ptr(), c_str(name), value )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_real32

function atlas_Config__get_real64(this, name, value) result(found)
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(inout) :: value
  integer :: found_int
  found_int = atlas__Config__get_double(this%c_ptr(), c_str(name), value )
  found = .False.
  if (found_int == 1) found = .True.
end function atlas_Config__get_real64

function atlas_Config__get_string(this, name, value) result(found)
  use atlas_Config_c_binding
  logical :: found
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(inout) :: value
  type(c_ptr) :: value_cptr
  integer :: found_int
  integer(c_int) :: value_size
  integer(c_int) :: value_allocated
  found_int = atlas__Config__get_string(this%c_ptr(),c_str(name),value_cptr,value_size,value_allocated)
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
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_int), intent(in) :: value(:)
  call atlas__Config__set_array_int(this%c_ptr(), c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_int32

subroutine atlas_Config__set_array_int64(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  integer(c_long), intent(in) :: value(:)
  call atlas__Config__set_array_long(this%c_ptr(), c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_int64

subroutine atlas_Config__set_array_real32(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value(:)
  call atlas__Config__set_array_float(this%c_ptr(), c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_real32

subroutine atlas_Config__set_array_real64(this, name, value)
  use atlas_Config_c_binding
  class(atlas_Config), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value(:)
  call atlas__Config__set_array_double(this%c_ptr(), c_str(name), &
    & value, size(value) )
end subroutine atlas_Config__set_array_real64

function atlas_Config__get_array_int32(this, name, value) result(found)
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
  found_int = atlas__Config__get_array_int(this%c_ptr(), c_str(name), &
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
  found_int = atlas__Config__get_array_long(this%c_ptr(), c_str(name), &
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
  found_int = atlas__Config__get_array_float(this%c_ptr(), c_str(name), &
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
  found_int = atlas__Config__get_array_double(this%c_ptr(), c_str(name), &
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
  use atlas_Config_c_binding
  character(len=:), allocatable :: json
  class(atlas_Config), intent(in) :: this
  type(c_ptr) :: json_cptr
  integer(c_int) :: json_size
  integer(c_int) :: json_allocated
  call atlas__Config__json(this%c_ptr(),json_cptr,json_size,json_allocated)
  allocate(character(len=json_size) :: json )
  json = c_ptr_to_string(json_cptr)
  if( json_allocated == 1 ) call c_ptr_free(json_cptr)
end function atlas_Config__json

end module atlas_config_module

