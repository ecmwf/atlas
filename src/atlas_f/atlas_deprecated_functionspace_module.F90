
module atlas_deprecated_functionspace_module

use iso_c_binding, only : c_ptr, c_int, c_long, c_double, c_float, c_f_pointer
use atlas_c_interop, only: c_str, c_to_f_string_cptr, view1d
use atlas_object_module, only : atlas_object
use atlas_metadata_module, only : atlas_Metadata
use atlas_field_module, only : atlas_Field
use atlas_gatherscatter_module, only: atlas_GatherScatter
use atlas_checksum_module, only: atlas_Checksum
use atlas_haloexchange_module, only: atlas_HaloExchange
implicit none

private :: c_ptr, c_int, c_long, c_double, c_float, c_f_pointer
private :: c_str, c_to_f_string_cptr, view1d
private :: atlas_object
private :: atlas_Metadata
private :: atlas_Field
private :: atlas_GatherScatter
private :: atlas_Checksum
private :: atlas_HaloExchange

public :: atlas_deprecated_FunctionSpace

private


!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_deprecated_FunctionSpace

! Purpose :
! -------
!   *FunctionSpace* :
!       Container type of fields that are defined on the same points
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done
!       Describes interpolation between nodes

! Methods :
! -------
!   name : The name or tag this function space was created with
!   create_field : Create a new real field in this function space with given name
!   remove_field : Remove a field with given name
!   field : Access to a field with given name
!   parallelise : Setup halo-exchange information
!   halo_exchange : Perform halo exchange on field_data

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: name => FunctionSpace__name
  procedure :: dof => FunctionSpace__dof
  procedure :: glb_dof => FunctionSpace__glb_dof
  procedure :: create_field => FunctionSpace__create_field
  procedure :: remove_field => FunctionSpace__remove_field
  procedure :: shape_arr => FunctionSpace__shape_arr
  procedure :: shape_idx => FunctionSpace__shape_idx
  generic :: shape => shape_arr, shape_idx
  procedure :: field => FunctionSpace__field
  procedure, public :: has_field => FunctionSpace__has_field
  procedure :: metadata => FunctionSpace__metadata
  procedure :: parallelise => FunctionSpace__parallelise
  procedure, private :: FunctionSpace__halo_exchange_int32_r1
  procedure, private :: FunctionSpace__halo_exchange_int32_r2
  procedure, private :: FunctionSpace__halo_exchange_int32_r3
  procedure, private :: FunctionSpace__halo_exchange_real32_r1
  procedure, private :: FunctionSpace__halo_exchange_real32_r2
  procedure, private :: FunctionSpace__halo_exchange_real32_r3
  procedure, private :: FunctionSpace__halo_exchange_real64_r1
  procedure, private :: FunctionSpace__halo_exchange_real64_r2
  procedure, private :: FunctionSpace__halo_exchange_real64_r3
  procedure, private :: FunctionSpace__halo_exchange_real64_r4
  procedure :: get_halo_exchange => FunctionSpace__get_halo_exchange
  procedure :: get_gather => FunctionSpace__get_gather
  procedure :: get_checksum => FunctionSpace__get_checksum
  generic :: halo_exchange => &
      & FunctionSpace__halo_exchange_int32_r1, &
      & FunctionSpace__halo_exchange_int32_r2, &
      & FunctionSpace__halo_exchange_int32_r3, &
      & FunctionSpace__halo_exchange_real32_r1, &
      & FunctionSpace__halo_exchange_real32_r2, &
      & FunctionSpace__halo_exchange_real32_r3, &
      & FunctionSpace__halo_exchange_real64_r1, &
      & FunctionSpace__halo_exchange_real64_r2, &
      & FunctionSpace__halo_exchange_real64_r3, &
      & FunctionSpace__halo_exchange_real64_r4
      procedure, private :: FunctionSpace__gather_real32_r1
      procedure, private :: FunctionSpace__gather_real32_r2
      procedure, private :: FunctionSpace__gather_real32_r3
  procedure, private :: FunctionSpace__gather_real64_r1
  procedure, private :: FunctionSpace__gather_real64_r2
  procedure, private :: FunctionSpace__gather_real64_r3
  procedure, private :: FunctionSpace__gather_int32_r1
  procedure, private :: FunctionSpace__gather_int32_r2
  generic :: gather => &
      & FunctionSpace__gather_real32_r1, &
      & FunctionSpace__gather_real32_r2, &
      & FunctionSpace__gather_real32_r3, &
      & FunctionSpace__gather_real64_r1, &
      & FunctionSpace__gather_real64_r2, &
      & FunctionSpace__gather_real64_r3, &
      & FunctionSpace__gather_int32_r1, &
      & FunctionSpace__gather_int32_r2

  procedure, public :: delete => atlas_deprecated_FunctionSpace__delete

END TYPE atlas_deprecated_FunctionSpace


!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

! -----------------------------------------------------------------------------
! FunctionSpace routines

function FunctionSpace__metadata(this) result(metadata)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  type(atlas_Metadata) :: Metadata
  call metadata%reset_c_ptr( atlas__deprecated__FunctionSpace__metadata(this%c_ptr()) )
end function FunctionSpace__metadata

subroutine FunctionSpace__create_field(this,name,nvars,kind)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: nvars
  integer, intent(in), optional :: kind

  ! ----------------------------------------------------
  ! ENUM DataType
  integer, parameter :: ATLAS_KIND_INT32  = -4
  integer, parameter :: ATLAS_KIND_INT64  = -8
  integer, parameter :: ATLAS_KIND_REAL32 =  4
  integer, parameter :: ATLAS_KIND_REAL64 =  8
  ! ----------------------------------------------------


  if (present(kind)) then
    if (kind == ATLAS_KIND_REAL64) then
      call atlas__deprecated__FunctionSpace__create_field_double(this%c_ptr(),c_str(name),nvars)
    else if (kind == ATLAS_KIND_REAL32) then
      call atlas__deprecated__FunctionSpace__create_field_float(this%c_ptr(),c_str(name),nvars)
    else if (kind == ATLAS_KIND_INT32) then
      call atlas__deprecated__FunctionSpace__create_field_int(this%c_ptr(),c_str(name),nvars)
    else if (kind == ATLAS_KIND_INT64) then
      call atlas__deprecated__FunctionSpace__create_field_long(this%c_ptr(),c_str(name),nvars)
    else
      write(0,*) "Unsupported kind"
      write(0,*) 'call abort()'
    endif
  else
    call atlas__deprecated__FunctionSpace__create_field_double(this%c_ptr(),c_str(name),nvars)
  end if
end subroutine FunctionSpace__create_field


subroutine FunctionSpace__remove_field(this,name)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__deprecated__FunctionSpace__remove_field(this%c_ptr(),c_str(name))
end subroutine FunctionSpace__remove_field

function FunctionSpace__name(this) result(name)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__deprecated__FunctionSpace__name(this%c_ptr())
  name = c_to_f_string_cptr(name_c_str)
end function FunctionSpace__name

function FunctionSpace__dof(this) result(dof)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer :: dof
  dof = atlas__deprecated__FunctionSpace__dof(this%c_ptr())
end function FunctionSpace__dof

function FunctionSpace__glb_dof(this) result(glb_dof)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer :: glb_dof
  glb_dof = atlas__deprecated__FunctionSpace__glb_dof(this%c_ptr())
end function FunctionSpace__glb_dof

function FunctionSpace__shape_arr(this) result(shape)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, pointer :: shape(:)
  type(c_ptr) :: shape_c_ptr
  integer(c_int) :: field_rank
  call atlas__deprecated__FunctionSpace__shapef(this%c_ptr(), shape_c_ptr, field_rank)
  call C_F_POINTER ( shape_c_ptr , shape , (/field_rank/) )
end function FunctionSpace__shape_arr

function FunctionSpace__shape_idx(this,idx) result(shape_idx)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, intent(in) :: idx
  integer :: shape_idx
  integer, pointer :: shape(:)
  shape => this%shape_arr()
  shape_idx = shape(idx)
end function FunctionSpace__shape_idx

function FunctionSpace__field(this,name) result(field)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: field
  field = atlas_Field( atlas__deprecated__FunctionSpace__field(this%c_ptr(), c_str(name) ) )
  if( field%is_null() ) write(0,*) 'call abort()'
  call field%return()
end function FunctionSpace__field

function FunctionSpace__has_field(this,name) result(flag)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__deprecated__FunctionSpace__has_field(this%c_ptr(), c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FunctionSpace__has_field

subroutine FunctionSpace__parallelise(this)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  call atlas__deprecated__FunctionSpace__parallelise(this%c_ptr())
end subroutine FunctionSpace__parallelise

function FunctionSpace__get_halo_exchange(this) result(halo_exchange)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  type(atlas_HaloExchange) :: halo_exchange
  call halo_exchange%reset_c_ptr( atlas__deprecated__FunctionSpace__halo_exchange( this%c_ptr() ) )
end function FunctionSpace__get_halo_exchange

subroutine FunctionSpace__halo_exchange_int32_r1(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__deprecated__FunctionSpace__halo_exchange_int( this%c_ptr(), field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r1
subroutine FunctionSpace__halo_exchange_int32_r2(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_int( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r2
subroutine FunctionSpace__halo_exchange_int32_r3(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_int( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r3

subroutine FunctionSpace__halo_exchange_real32_r1(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__deprecated__FunctionSpace__halo_exchange_float( this%c_ptr(), field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r1
subroutine FunctionSpace__halo_exchange_real32_r2(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_float( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r2
subroutine FunctionSpace__halo_exchange_real32_r3(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_float( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r3

subroutine FunctionSpace__halo_exchange_real64_r1(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__deprecated__FunctionSpace__halo_exchange_double( this%c_ptr(), field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r1
subroutine FunctionSpace__halo_exchange_real64_r2(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_double( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r2
subroutine FunctionSpace__halo_exchange_real64_r3(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_double( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r3
subroutine FunctionSpace__halo_exchange_real64_r4(this, field_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__deprecated__FunctionSpace__halo_exchange_double( this%c_ptr(), view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r4


function FunctionSpace__get_gather(this) result(gather)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  type(atlas_GatherScatter) :: gather
  call gather%reset_c_ptr( atlas__deprecated__FunctionSpace__gather( this%c_ptr() ) )
end function FunctionSpace__get_gather

subroutine FunctionSpace__gather_real32_r1(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_float), intent(in) :: field_data(:)
  real(c_float), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__deprecated__FunctionSpace__gather_float( this%c_ptr(), field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r1


subroutine FunctionSpace__gather_real32_r2(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__deprecated__FunctionSpace__gather_float( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r2


subroutine FunctionSpace__gather_real32_r3(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__deprecated__FunctionSpace__gather_float( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r3

subroutine FunctionSpace__gather_real64_r1(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(in) :: field_data(:)
  real(c_double), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__deprecated__FunctionSpace__gather_double( this%c_ptr(), field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r1


subroutine FunctionSpace__gather_real64_r2(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__deprecated__FunctionSpace__gather_double( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r2


subroutine FunctionSpace__gather_real64_r3(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__deprecated__FunctionSpace__gather_double( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r3

subroutine FunctionSpace__gather_int32_r1(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, intent(in) :: field_data(:)
  integer, intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__deprecated__FunctionSpace__gather_int( this%c_ptr(), field_data, size(field_data), &
                                       & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r1

subroutine FunctionSpace__gather_int32_r2(this, field_data, glbfield_data)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  integer, intent(in) :: field_data(:,:)
  integer, intent(inout) :: glbfield_data(:,:)
  integer, pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__deprecated__FunctionSpace__gather_int( this%c_ptr(), view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r2

function FunctionSpace__get_checksum(this) result(checksum)
  use atlas_FunctionSpace_c_binding
  class(atlas_deprecated_FunctionSpace), intent(in) :: this
  type(atlas_Checksum) :: checksum
  call checksum%reset_c_ptr( atlas__deprecated__FunctionSpace__checksum( this%c_ptr() ) )
end function FunctionSpace__get_checksum

subroutine atlas_deprecated_FunctionSpace__delete(this)
  class(atlas_deprecated_FunctionSpace), intent(inout) :: this
end subroutine

end module atlas_deprecated_functionspace_module

