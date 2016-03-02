
module atlas_Value_module

use, intrinsic :: iso_c_binding, only : c_int, c_long, c_double, c_float, c_ptr
use atlas_c_interop, only : c_str, c_to_f_string_cptr
use atlas_object_module, only : atlas_object
implicit none

private :: c_int, c_long, c_double, c_float, c_ptr
private :: c_str, c_to_f_string_cptr

public :: atlas_Value

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_Value

! Purpose :
! -------
!   *Value* : Container of Value, which internally stores any standard type

! Methods :
! -------
!   get : Return a value in standard type

! Author :
! ------
!   June-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: get_logical => atlas_Value__get_logical
  procedure :: get_int32   => atlas_Value__get_int32
  procedure :: get_int64   => atlas_Value__get_int64
  procedure :: get_real32  => atlas_Value__get_real32
  procedure :: get_real64  => atlas_Value__get_real64
  procedure :: get_string  => atlas_Value__get_string
  procedure :: get_array_int32  => atlas_Value__get_array_int32
  procedure :: get_array_int64  => atlas_Value__get_array_int64
  procedure :: get_array_real32 => atlas_Value__get_array_real32
  procedure :: get_array_real64 => atlas_Value__get_array_real64
  generic :: get => get_logical, get_int32, get_int64, get_real32, get_real64, get_string, &
                    get_array_int32, get_array_int64, get_array_real32, get_array_real64


  procedure, public :: delete => atlas_Value__delete

END TYPE atlas_Value

!------------------------------------------------------------------------------

interface atlas_Value
  module procedure atlas_Value__ctor_int32
  module procedure atlas_Value__ctor_int64
  module procedure atlas_Value__ctor_real32
  module procedure atlas_Value__ctor_real64
  module procedure atlas_Value__ctor_string
  module procedure atlas_Value__ctor_array_int32
  module procedure atlas_Value__ctor_array_int64
  module procedure atlas_Value__ctor_array_real32
  module procedure atlas_Value__ctor_array_real64
end interface

!------------------------------------------------------------------------------


! =============================================================================
CONTAINS
! =============================================================================

! -----------------------------------------------------------------------------
! Value routines

function atlas_Value__ctor_int32(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_int
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  integer(c_int), intent(in) :: val
  call Value%reset_c_ptr( atlas__Value__new_int(val) )
end function

function atlas_Value__ctor_int64(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_long
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  integer(c_long), intent(in) :: val
  call Value%reset_c_ptr( atlas__Value__new_long(val) )
end function

function atlas_Value__ctor_real32(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_float
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  real(c_float), intent(in) :: val
  call Value%reset_c_ptr( atlas__Value__new_float(val) )
end function

function atlas_Value__ctor_real64(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_double
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  real(c_double), intent(in) :: val
  call Value%reset_c_ptr( atlas__Value__new_double(val) )
end function

function atlas_Value__ctor_string(val) result(Value)
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  character(len=*), intent(in) :: val
  call Value%reset_c_ptr( atlas__Value__new_string(c_str(val)) )
end function

function atlas_Value__ctor_array_int32(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_int
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  integer(c_int), intent(in) :: val(:)
  call Value%reset_c_ptr( atlas__Value__new_array_int(val,size(val)) )
end function

function atlas_Value__ctor_array_int64(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_long
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  integer(c_long), intent(in) :: val(:)
  call Value%reset_c_ptr( atlas__Value__new_array_long(val,size(val)) )
end function

function atlas_Value__ctor_array_real32(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_float
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  real(c_float), intent(in) :: val(:)
  call Value%reset_c_ptr( atlas__Value__new_array_float(val,size(val)) )
end function

function atlas_Value__ctor_array_real64(val) result(Value)
  use, intrinsic :: iso_c_binding , only : c_double
  use atlas_atlas_value_c_binding
  type(atlas_Value) :: Value
  real(c_double), intent(in) :: val(:)
  call Value%reset_c_ptr( atlas__Value__new_array_double(val,size(val)) )
end function

subroutine atlas_Value__delete(this)
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Value__delete(this%c_ptr())
  endif
  call this%reset_c_ptr()
end subroutine atlas_Value__delete


subroutine atlas_Value__array_delete(this)
  use atlas_atlas_value_c_binding
  type(atlas_Value), intent(inout) :: this(:)
  integer :: i
  do i=lbound(this,1),ubound(this,1)
    call atlas_Value__delete(this(i))
  enddo
end subroutine atlas_Value__array_delete



subroutine atlas_Value__get_logical(this,val)
  use, intrinsic :: iso_c_binding , only : c_int
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  logical, intent(out) :: val
  integer(c_int) :: ival
  call atlas__Value__int(this%c_ptr(),ival)
  if( ival == 0 ) then
    val = .False.
  else
    val = .True.
  endif
end subroutine

subroutine atlas_Value__get_int32(this,val)
  use, intrinsic :: iso_c_binding , only : c_int
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  integer(c_int), intent(out) :: val
  call atlas__Value__int(this%c_ptr(),val)
end subroutine

subroutine atlas_Value__get_int64(this,val)
  use, intrinsic :: iso_c_binding , only : c_long
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  integer(c_long), intent(out) :: val
  call atlas__Value__long(this%c_ptr(),val)
end subroutine

subroutine atlas_Value__get_real32(this,val)
  use, intrinsic :: iso_c_binding , only : c_float
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  real(c_float), intent(out) :: val
  call atlas__Value__float(this%c_ptr(),val)
end subroutine

subroutine atlas_Value__get_real64(this,val)
  use, intrinsic :: iso_c_binding , only : c_double
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  real(c_double), intent(out) :: val
  call atlas__Value__double(this%c_ptr(),val)
end subroutine

subroutine atlas_Value__get_string(this,val)
  use, intrinsic :: iso_c_binding , only : c_ptr, c_int
  use atlas_c_interop
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  character(len=:), allocatable, intent(out) :: val
  type(c_ptr) :: val_cptr
  integer(c_int) :: val_size
  integer(c_int) :: val_allocated
  call atlas__Value__string(this%c_ptr(),val_cptr,val_size,val_allocated)
  val = c_to_f_string_cptr(val_cptr)
  if( val_allocated == 1 ) call atlas_free(val_cptr)
end subroutine

subroutine atlas_Value__get_array_int32(this,val)
  use, intrinsic :: iso_c_binding , only : c_int, c_f_pointer, c_ptr
  use atlas_c_interop
  use atlas_atlas_value_c_binding
  class(atlas_Value), intent(in) :: this
  integer(c_int), allocatable, intent(out) :: val(:)
  integer(c_int), pointer :: val_fptr(:)
  type(c_ptr) :: val_cptr
  integer(c_int) :: val_size
  integer(c_int) :: val_allocated
  call atlas__Value__array_int(this%c_ptr(),val_cptr,val_size,val_allocated)
  call c_f_pointer(val_cptr,val_fptr,(/val_size/))
  if( allocated(val) ) deallocate(val)
  allocate( val(val_size) )
  val(:) = val_fptr(:)
  if( val_allocated == 1 ) call atlas_free(val_cptr)
end subroutine

subroutine atlas_Value__get_array_int64(this,val)
  use, intrinsic :: iso_c_binding , only : c_long, c_int, c_ptr, c_f_pointer
  use atlas_atlas_value_c_binding
  use atlas_c_interop
  class(atlas_Value), intent(in) :: this
  integer(c_long), allocatable, intent(out) :: val(:)
  integer(c_long), pointer :: val_fptr(:)
  type(c_ptr) :: val_cptr
  integer(c_int) :: val_size
  integer(c_int) :: val_allocated
  call atlas__Value__array_long(this%c_ptr(),val_cptr,val_size,val_allocated)
  call c_f_pointer(val_cptr,val_fptr,(/val_size/))
  if( allocated(val) ) deallocate(val)
  allocate( val(val_size) )
  val(:) = val_fptr(:)
  if( val_allocated == 1 ) call atlas_free(val_cptr)
end subroutine

subroutine atlas_Value__get_array_real32(this,val)
  use, intrinsic :: iso_c_binding , only : c_float, c_int, c_ptr, c_f_pointer
  use atlas_atlas_value_c_binding
  use atlas_c_interop
  class(atlas_Value), intent(in) :: this
  real(c_float), allocatable, intent(out) :: val(:)
  real(c_float), pointer :: val_fptr(:)
  type(c_ptr) :: val_cptr
  integer(c_int) :: val_size
  integer(c_int) :: val_allocated
  call atlas__Value__array_float(this%c_ptr(),val_cptr,val_size,val_allocated)
  call c_f_pointer(val_cptr,val_fptr,(/val_size/))
  if( allocated(val) ) deallocate(val)
  allocate( val(val_size) )
  val(:) = val_fptr(:)
  if( val_allocated == 1 ) call atlas_free(val_cptr)
end subroutine

subroutine atlas_Value__get_array_real64(this,val)
  use, intrinsic :: iso_c_binding , only : c_double, c_int, c_ptr, c_f_pointer
  use atlas_atlas_value_c_binding
  use atlas_c_interop
  class(atlas_Value), intent(in) :: this
  real(c_double), allocatable, intent(out) :: val(:)
  real(c_double), pointer :: val_fptr(:)
  type(c_ptr) :: val_cptr
  integer(c_int) :: val_size
  integer(c_int) :: val_allocated
  call atlas__Value__array_double(this%c_ptr(),val_cptr,val_size,val_allocated)
  call c_f_pointer(val_cptr,val_fptr,(/val_size/))
  if( allocated(val) ) deallocate(val)
  allocate( val(val_size) )
  val(:) = val_fptr(:)
  if( val_allocated == 1 ) call atlas_free(val_cptr)
end subroutine

end module atlas_Value_module

