! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_allocate_module

implicit none
private

public :: atlas_allocate_managedmem
public :: atlas_deallocate_managedmem


interface atlas_allocate_managedmem
  module procedure atlas_allocate_managedmem_real64_r1
  module procedure atlas_allocate_managedmem_real32_r1
  module procedure atlas_allocate_managedmem_int32_r1
  module procedure atlas_allocate_managedmem_int64_r1_int32
  module procedure atlas_allocate_managedmem_int64_r1_int64
  module procedure atlas_allocate_managedmem_real64_r2
  module procedure atlas_allocate_managedmem_real32_r2
  module procedure atlas_allocate_managedmem_int32_r2
  module procedure atlas_allocate_managedmem_int64_r2_int32
  module procedure atlas_allocate_managedmem_int64_r2_int64
end interface

interface atlas_deallocate_managedmem
  module procedure atlas_deallocate_managedmem_real64_r1
  module procedure atlas_deallocate_managedmem_real32_r1
  module procedure atlas_deallocate_managedmem_int32_r1
  module procedure atlas_deallocate_managedmem_int64_r1
  module procedure atlas_deallocate_managedmem_real64_r2
  module procedure atlas_deallocate_managedmem_real32_r2
  module procedure atlas_deallocate_managedmem_int32_r2
  module procedure atlas_deallocate_managedmem_int64_r2
end interface

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

subroutine atlas_allocate_managedmem_real64_r1( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_double), pointer :: a(:)
  integer(c_int) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)
  if( size > 0 ) then
    call atlas__allocate_managedmem_double( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_real32_r1( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_float), pointer :: a(:)
  integer(c_int) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)
  if( size > 0 ) then
    call atlas__allocate_managedmem_float( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_int32_r1( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_int), pointer :: a(:)
  integer(c_int) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)
  if( size > 0 ) then
    call atlas__allocate_managedmem_int( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_int64_r1_int32( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_long), pointer :: a(:)
  integer(c_int32_t) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)
  if( size > 0 ) then
    call atlas__allocate_managedmem_long( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_int64_r1_int64( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_long), pointer :: a(:)
  integer(c_int64_t) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)
  if( size > 0 ) then
    call atlas__allocate_managedmem_long( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_int32_r2( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_int), pointer :: a(:,:)
  integer(c_int) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)*dims(2)
  if( size > 0 ) then
    call atlas__allocate_managedmem_int( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_int64_r2_int32( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_long), pointer :: a(:,:)
  integer(c_int32_t) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)*dims(2)
  if( size > 0 ) then
    call atlas__allocate_managedmem_long( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_int64_r2_int64( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_long), pointer :: a(:,:)
  integer(c_int64_t) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)*dims(2)
  if( size > 0 ) then
    call atlas__allocate_managedmem_long( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_real64_r2( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_double), pointer :: a(:,:)
  integer(c_int) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)*dims(2)
  if( size > 0 ) then
    call atlas__allocate_managedmem_double( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

subroutine atlas_allocate_managedmem_real32_r2( A, dims )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_float), pointer :: a(:,:)
  integer(c_int) :: dims(:)
  type(c_ptr) :: value_cptr
  integer(c_size_t) :: size
  size = dims(1)*dims(2)
  if( size > 0 ) then
    call atlas__allocate_managedmem_float( value_cptr, size )
    call c_f_pointer(value_cptr,a,dims)
  endif
end subroutine

!------------------------------------------------------------------------------

! These functions are private in fckit_array_module

function c_loc_int32(x)
  use, intrinsic :: iso_c_binding
  integer(c_int32_t), target :: x
  type(c_ptr) :: c_loc_int32
  c_loc_int32 = c_loc(x)
end function

! =============================================================================

function c_loc_int64(x)
  use, intrinsic :: iso_c_binding
  integer(c_int64_t), target :: x
  type(c_ptr) :: c_loc_int64
  c_loc_int64 = c_loc(x)
end function

! =============================================================================

function c_loc_real32(x)
  use, intrinsic :: iso_c_binding
  real(c_float), target :: x
  type(c_ptr) :: c_loc_real32
  c_loc_real32 = c_loc(x)
end function

! =============================================================================

function c_loc_real64(x)
  use, intrinsic :: iso_c_binding
  real(c_double), target :: x
  type(c_ptr) :: c_loc_real64
  c_loc_real64 = c_loc(x)
end function

! =============================================================================

!------------------------------------------------------------------------------

subroutine atlas_deallocate_managedmem_real64_r1( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_double), pointer :: a(:)
  call atlas__deallocate_managedmem_double( c_loc_real64(A(1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_real32_r1( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_float), pointer :: a(:)
  call atlas__deallocate_managedmem_float( c_loc_real32(A(1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_int32_r1( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_int), pointer :: a(:)
  call atlas__deallocate_managedmem_int( c_loc_int32(A(1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_int64_r1( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_long), pointer :: a(:)
  call atlas__deallocate_managedmem_long( c_loc_int64(A(1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_real64_r2( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_double), pointer :: a(:,:)
  call atlas__deallocate_managedmem_double( c_loc_real64(A(1,1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_real32_r2( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  real(c_float), pointer :: a(:,:)
  call atlas__deallocate_managedmem_float( c_loc_real32(A(1,1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_int32_r2( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_int), pointer :: a(:,:)
  call atlas__deallocate_managedmem_int( c_loc_int32(A(1,1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

subroutine atlas_deallocate_managedmem_int64_r2( A )
  use, intrinsic :: iso_c_binding
  use atlas_allocate_c_binding
  integer(c_long), pointer :: a(:,:)
  call atlas__deallocate_managedmem_long( c_loc_int64(A(1,1)), size(A,KIND=c_size_t) )
  nullify( a )
end subroutine

!------------------------------------------------------------------------------

end module atlas_allocate_module

