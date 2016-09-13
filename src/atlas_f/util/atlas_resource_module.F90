
module atlas_resource_module

use, intrinsic :: iso_c_binding, only : c_char, c_int, c_long, c_double, c_float, c_ptr
use fckit_c_interop_module, only : c_str, c_ptr_to_string
implicit none

private :: c_char, c_int, c_long, c_double, c_float, c_ptr
private :: c_str, c_ptr_to_string

public :: atlas_resource
public :: atlas_resource_set

private


!------------------------------------------------------------------------------

INTERFACE atlas_resource

! Purpose :
! -------
!   *resource* : Configuration
!
! Author :
! ------
!   20-dec-2014 Willem Deconinck     *ECMWF*
! -----------------------------------------------------------------------------
  module procedure resource_get_int32
  module procedure resource_get_int64
  module procedure resource_get_real32
  module procedure resource_get_real64
  module procedure resource_get_string
end interface atlas_resource

!------------------------------------------------------------------------------

INTERFACE atlas_resource_set

! Purpose :
! -------
!   *resource* : Configuration
!
! Author :
! ------
!   10-june-2015 Willem Deconinck     *ECMWF*
! -----------------------------------------------------------------------------
  module procedure resource_set_int32
  module procedure resource_set_int64
  module procedure resource_set_real32
  module procedure resource_set_real64
  module procedure resource_set_string
end interface atlas_resource_set

! =============================================================================
CONTAINS
! =============================================================================


subroutine resource_get_int32(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_int), intent(in) :: default_value
  integer(c_int), intent(out) :: value
  value = atlas__resource_int( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_int64(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_long), intent(in) :: default_value
  integer(c_long), intent(out) :: value
  value = atlas__resource_long( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_real32(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_float), intent(in) :: default_value
  real(c_float), intent(out) :: value
  value = atlas__resource_float( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_real64(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_double), intent(in) :: default_value
  real(c_double), intent(out) :: value
  value = atlas__resource_double( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_string(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  character(len=*), intent(in) :: default_value
  character(len=*), intent(out) :: value
  type(c_ptr) :: value_c_str
  value_c_str = atlas__resource_string( c_str(resource_str), c_str(default_value) )
  value = c_ptr_to_string(value_c_str)
end subroutine

subroutine resource_set_int32(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_int), intent(in) :: value
  call atlas__resource_set_int( c_str(resource_str), value )
end subroutine

subroutine resource_set_int64(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_long), intent(in) ::value
  call atlas__resource_set_long( c_str(resource_str), value )
end subroutine

subroutine resource_set_real32(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_float), intent(in) :: value
  call atlas__resource_set_float( c_str(resource_str), value )
end subroutine

subroutine resource_set_real64(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_double), intent(in) :: value
  call atlas__resource_set_double( c_str(resource_str), value )
end subroutine

subroutine resource_set_string(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  character(len=*), intent(in) :: value
  type(c_ptr) :: value_c_str
  call atlas__resource_set_string( c_str(resource_str), c_str(value) )
end subroutine


end module atlas_resource_module

