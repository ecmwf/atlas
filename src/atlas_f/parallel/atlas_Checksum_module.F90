#include "atlas/atlas_f.h"

module atlas_checksum_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_long, c_float, c_double, c_char
use fckit_array_module, only : array_stride, array_view1d
use fckit_c_interop_module, only : c_str_to_string
use fckit_object_module, only : fckit_object

implicit none

private :: c_ptr, c_int, c_long, c_float, c_double, c_char
private :: array_stride, array_view1d, c_str_to_string
private :: fckit_object

public :: atlas_Checksum

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_object) :: atlas_Checksum

! Purpose :
! -------
!   *Checksum* :

! Methods :
! -------
!   setup : Setup using arrays detailing proc, glb_idx, remote_idx, max_glb_idx
!   execute : Do the Checksum

! Author :
! ------
!   27-Jun-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: Checksum__setup32
  procedure, private :: Checksum__setup64
  procedure, private :: Checksum__execute_int32_r1
  procedure, private :: Checksum__execute_int32_r2
  procedure, private :: Checksum__execute_int32_r3
  procedure, private :: Checksum__execute_real32_r1
  procedure, private :: Checksum__execute_real32_r2
  procedure, private :: Checksum__execute_real32_r3
  procedure, private :: Checksum__execute_real64_r1
  procedure, private :: Checksum__execute_real64_r3
  procedure, private :: Checksum__execute_real64_r2
  generic :: setup => &
      & Checksum__setup32, &
      & Checksum__setup64
  generic :: execute => &
      & Checksum__execute_int32_r1, &
      & Checksum__execute_int32_r2, &
      & Checksum__execute_int32_r3, &
      & Checksum__execute_real32_r1, &
      & Checksum__execute_real32_r2, &
      & Checksum__execute_real32_r3, &
      & Checksum__execute_real64_r1, &
      & Checksum__execute_real64_r2, &
      & Checksum__execute_real64_r3

  procedure, public :: delete => atlas_Checksum__delete

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Checksum__final_auto
#endif

END TYPE atlas_Checksum

!------------------------------------------------------------------------------

interface atlas_Checksum
  module procedure atlas_Checksum__ctor
end interface

!------------------------------------------------------------------------------
!========================================================
contains
!========================================================

! ------------------------------------------------------------------------------
! Checksum routines

function atlas_Checksum__ctor() result(Checksum)
  use atlas_checksum_c_binding
  type(atlas_Checksum) :: Checksum
  call Checksum%reset_c_ptr( atlas__Checksum__new() )
end function atlas_checksum__ctor

subroutine atlas_Checksum__delete(this)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Checksum__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_Checksum__delete

subroutine Checksum__setup32(this, part, remote_idx, glb_idx )
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  integer(c_int), intent(in) :: part(:)
  integer(c_int), intent(in) :: remote_idx(:)
  integer(c_int), intent(in) :: glb_idx(:)
  call atlas__Checksum__setup32( this%c_ptr(), part, remote_idx, 1, &
    &                          glb_idx, size(part) )
end subroutine Checksum__setup32

subroutine Checksum__setup64(this, part, remote_idx, glb_idx)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  integer(c_int), intent(in) :: part(:)
  integer(c_int), intent(in) :: remote_idx(:)
  integer(c_long), intent(in) :: glb_idx(:)
  call atlas__Checksum__setup64( this%c_ptr(), part, remote_idx, 1, &
    &                            glb_idx, size(part) )
end subroutine Checksum__setup64

function Checksum__execute_int32_r1(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer :: lstrides(1), lextents(1), lrank=1
  lstrides = (/ array_stride(loc_field_data,2) /)
  lextents = (/ 1                        /)
  call atlas__Checksum__execute_strided_int( this%c_ptr(), &
    &  loc_field_data, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_int32_r1

function Checksum__execute_int32_r2(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer, pointer :: lview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  lstrides = (/ array_stride(loc_field_data,2), array_stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_int( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_int32_r2


function Checksum__execute_int32_r3(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  integer, intent(in)  :: loc_field_data(:,:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer, pointer :: lview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  lstrides = (/ array_stride(loc_field_data,3), array_stride(loc_field_data,2) , array_stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size(loc_field_data,1) /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_int( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_int32_r3

function Checksum__execute_real32_r1(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  real(c_float), intent(in)   :: loc_field_data(:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer :: lstrides(1), lextents(1), lrank=1
  lstrides = (/ array_stride(loc_field_data,1) /)
  lextents = (/ 1                        /)
  call atlas__Checksum__execute_strided_float( this%c_ptr(), &
    &  loc_field_data, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_real32_r1
function Checksum__execute_real32_r2(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_float), pointer :: lview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  lstrides = (/ array_stride(loc_field_data,2), array_stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_float( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_real32_r2
function Checksum__execute_real32_r3(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  real(c_float), intent(in)  :: loc_field_data(:,:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_float), pointer :: lview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  lstrides = (/ array_stride(loc_field_data,3), array_stride(loc_field_data,2) , array_stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size  (loc_field_data,1) /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_float( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_real32_r3

function Checksum__execute_real64_r1(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  real(c_double), intent(in)   :: loc_field_data(:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  integer :: lstrides(1), lextents(1), lrank=1
  real(c_double), pointer :: lview(:)
  lstrides = (/ array_stride(loc_field_data,1) /)
  lextents = (/ 1                        /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_double( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_real64_r1
function Checksum__execute_real64_r2(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_double), pointer :: lview(:)
  integer :: lstrides(2), lextents(2), lrank=2
  lstrides = (/ array_stride(loc_field_data,2), array_stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,1) /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_double( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_real64_r2
function Checksum__execute_real64_r3(this, loc_field_data) result(checksum)
  use atlas_checksum_c_binding
  class(atlas_Checksum), intent(in) :: this
  real(c_double), intent(in)  :: loc_field_data(:,:,:)
  character(len=:), allocatable :: checksum
  character(kind=c_char) :: checksum_c_str(132)
  real(c_double), pointer :: lview(:)
  integer :: lstrides(3), lextents(3), lrank=3
  lstrides = (/ array_stride(loc_field_data,3), array_stride(loc_field_data,2) , array_stride(loc_field_data,1) /)
  lextents = (/ 1,                        size  (loc_field_data,2) , size  (loc_field_data,1) /)
  lview => array_view1d(loc_field_data)
  call atlas__Checksum__execute_strided_double( this%c_ptr(), &
    &  lview, lstrides, lextents, lrank, checksum_c_str )
  checksum = c_str_to_string(checksum_c_str)
end function Checksum__execute_real64_r3

!-------------------------------------------------------------------------------

subroutine atlas_Checksum__final_auto(this)
  type(atlas_Checksum) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Checksum__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! -----------------------------------------------------------------------------

end module atlas_checksum_module

