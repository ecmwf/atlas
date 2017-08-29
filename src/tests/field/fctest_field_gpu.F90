! (C) Copyright 1996-2016 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_Field_gpu_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none

contains

end module


subroutine test_res(ie, je, v1, vres)
  implicit none
  integer :: ie, je
  real(8), intent(in) :: v1(ie, je)
  real(8), intent(out) :: vres
  integer :: i,j

 !$acc kernels deviceptr(v1) copyout(vres)
     vres = v1(1,1)
 !$acc end kernels

end subroutine test_res


! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_Field_gpu,fcta_Field_gpu_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_host_data )
type(atlas_Field) :: field
real(8), pointer :: host(:,:)
real(8), pointer :: device_ptr(:,:)
real(8),pointer :: v1(:, :)
real(8) :: vres

field = atlas_Field(kind=atlas_real(8),shape=[10,5])

FCTEST_CHECK( .not. field%host_needs_update() )
FCTEST_CHECK( field%device_needs_update() )

call field%clone_to_device()
call field%host_data(host)
call field%device_data(device_ptr)

call test_res(10,5,device_ptr, vres)


!! acc kernels deviceptr(v1)  copyout(vres)
!  v1(2,2) = 3.5
!  vres = 2!v1(2,2)
!! acc end kernels

!FCTEST_CHECK_EQUAL( val, 3.5_c_float )


FCTEST_CHECK( .not. field%device_needs_update() )

call field%final()
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

