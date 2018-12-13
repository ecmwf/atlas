! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module allocate_unified_module

interface

  subroutine allocate_unified_impl( a, N ) bind(C,name="allocate_unified_impl")
    use, intrinsic :: iso_c_binding
    type(c_ptr) :: a
    integer(c_int), value :: N
  end subroutine

end interface

contains

  subroutine allocate_unified( a, N )
    use, intrinsic :: iso_c_binding
    real(c_double), pointer :: a(:)
    integer(c_int) :: N
    type(c_ptr) :: value_cptr
    call allocate_unified_impl(value_cptr, N)
    call c_f_pointer(value_cptr,a,(/N/))
  end subroutine

end module

TESTSUITE(fctest_atlas_openacc_with_unified_memory)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_openacc_with_unified_memory )
  use, intrinsic :: iso_c_binding
  use allocate_unified_module
  implicit none

  real(c_double), pointer :: value(:)
  integer :: i
  integer :: N
  N = 100

  call allocate_unified(value, N)

  ! On host:
  do i=1,N
    value(i) = i
  enddo

  ! On device (automatic copies on access because unified memory)
!$acc kernels
  do i=1,N
    value(i) = value(i) * 10
  enddo
!$acc end kernels

  ! On host:
  do i=1,N
    FCTEST_CHECK_CLOSE( value(i), real( 10*i, c_double), real( 0.0001, c_double ) )
  enddo
END_TEST

! -----------------------------------------------------------------------------
!
! NOTE: To be certain that cuda kernel is being generated and launched, run this
!       test with "nvprof -s <executable>" 
!
! -----------------------------------------------------------------------------

END_TESTSUITE

