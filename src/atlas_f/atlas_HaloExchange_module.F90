
module atlas_haloexchange_module

use, intrinsic :: iso_c_binding, only : c_ptr, c_long, c_float, c_double
use atlas_c_interop, only : stride, view1d
use atlas_object_module, only : atlas_object

implicit none

private :: c_ptr, c_long, c_float, c_double
private :: stride, view1d
private :: atlas_object

public :: atlas_HaloExchange

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_HaloExchange

! Purpose :
! -------
!   *HaloExchange* :

! Methods :
! -------
!   setup : Setup using arrays detailing proc and glb_idx, bounds and parbound
!   execute : Do the halo exchange

! Author :
! ------
!   17-Dec-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: setup => HaloExchange__setup
  procedure, private :: HaloExchange__execute_int32_r1
  procedure, private :: HaloExchange__execute_int32_r2
  procedure, private :: HaloExchange__execute_int32_r3
  procedure, private :: HaloExchange__execute_int64_r1
  procedure, private :: HaloExchange__execute_int64_r2
  procedure, private :: HaloExchange__execute_int64_r3
  procedure, private :: HaloExchange__execute_real32_r1
  procedure, private :: HaloExchange__execute_real32_r2
  procedure, private :: HaloExchange__execute_real32_r3
  procedure, private :: HaloExchange__execute_real64_r1
  procedure, private :: HaloExchange__execute_real64_r2
  procedure, private :: HaloExchange__execute_real64_r3
  procedure, private :: HaloExchange__execute_real64_r4
  generic :: execute => &
      & HaloExchange__execute_int32_r1, &
      & HaloExchange__execute_int32_r2, &
      & HaloExchange__execute_int32_r3, &
      & HaloExchange__execute_int64_r1, &
      & HaloExchange__execute_int64_r2, &
      & HaloExchange__execute_int64_r3, &
      & HaloExchange__execute_real32_r1, &
      & HaloExchange__execute_real32_r2, &
      & HaloExchange__execute_real32_r3, &
      & HaloExchange__execute_real64_r1, &
      & HaloExchange__execute_real64_r2, &
      & HaloExchange__execute_real64_r3, &
      & HaloExchange__execute_real64_r4

  procedure, public :: delete => atlas_HaloExchange__delete

END TYPE atlas_HaloExchange
!------------------------------------------------------------------------------

interface atlas_HaloExchange
  module procedure atlas_HaloExchange__ctor
end interface

!========================================================
contains
!========================================================


! ------------------------------------------------------------------------------
! HaloExchange routines

function atlas_HaloExchange__ctor() result(halo_exchange)
  use atlas_haloexchange_c_binding
  type(atlas_HaloExchange) :: halo_exchange
  call halo_exchange%reset_c_ptr( atlas__HaloExchange__new() )
end function atlas_HaloExchange__ctor

subroutine atlas_HaloExchange__delete(this)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__HaloExchange__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine atlas_HaloExchange__delete

subroutine HaloExchange__setup(this, part, remote_idx)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer, intent(in) :: part(:)
  integer, intent(in) :: remote_idx(:)
  call atlas__HaloExchange__setup( this%c_ptr(), part, remote_idx, 1, size(part) )
end subroutine HaloExchange__setup


subroutine HaloExchange__execute_int32_r1(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer, intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_int( this%c_ptr(), field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_int32_r1

subroutine HaloExchange__execute_int32_r2(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_int( this%c_ptr(), view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_int32_r2

subroutine HaloExchange__execute_int32_r3(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  integer :: strides(3), extents(3)
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_int( this%c_ptr(), view, &
    & strides, extents, 3 )
end subroutine HaloExchange__execute_int32_r3

subroutine HaloExchange__execute_int64_r1(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer(c_long), intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_long( this%c_ptr(), field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_int64_r1

subroutine HaloExchange__execute_int64_r2(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer(c_long), intent(inout) :: field_data(:,:)
  integer(c_long), pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_long( this%c_ptr(), view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_int64_r2

subroutine HaloExchange__execute_int64_r3(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  integer(c_long), intent(inout) :: field_data(:,:,:)
  integer(c_long), pointer :: view(:)
  integer :: strides(3), extents(3)
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_long( this%c_ptr(), view, &
    & strides, extents, 3 )
end subroutine HaloExchange__execute_int64_r3

subroutine HaloExchange__execute_real32_r1(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_float( this%c_ptr(), field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_real32_r1
subroutine HaloExchange__execute_real32_r2(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_float( this%c_ptr(), view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_real32_r2
subroutine HaloExchange__execute_real32_r3(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  integer :: strides(3), extents(3), rank=3
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_float( this%c_ptr(), view, &
    & strides, extents, rank )
end subroutine HaloExchange__execute_real32_r3

subroutine HaloExchange__execute_real64_r1(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
  integer :: strides(1), extents(1)
  strides = (/ stride(field_data,1) /)
  extents = (/ 1                    /)
  call atlas__HaloExchange__execute_strided_double( this%c_ptr(), field_data, &
    & strides, extents, 1 )
end subroutine HaloExchange__execute_real64_r1
subroutine HaloExchange__execute_real64_r2(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  integer :: strides(2), extents(2)
  view => view1d(field_data)
  strides = (/ stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1                    , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_double( this%c_ptr(), view, &
    & strides, extents, 2 )
end subroutine HaloExchange__execute_real64_r2
subroutine HaloExchange__execute_real64_r3(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  integer :: strides(3), extents(3), rank=3
  view => view1d(field_data)
  strides = (/ stride(field_data,3), stride(field_data,2) , stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,2) , ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_double( this%c_ptr(), view, &
    & strides, extents, rank )
end subroutine HaloExchange__execute_real64_r3
subroutine HaloExchange__execute_real64_r4(this, field_data)
  use atlas_haloexchange_c_binding
  class(atlas_HaloExchange), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  integer :: strides(4), extents(4), rank=4
  view => view1d(field_data)
  strides = (/ stride(field_data,4), stride(field_data,3), stride(field_data,2), stride(field_data,1) /)
  extents = (/ 1,                    ubound(field_data,3), ubound(field_data,2), ubound(field_data,1) /)
  call atlas__HaloExchange__execute_strided_double( this%c_ptr(), view, &
    & strides, extents, rank )
end subroutine HaloExchange__execute_real64_r4

! -----------------------------------------------------------------------------

end module atlas_haloexchange_module

