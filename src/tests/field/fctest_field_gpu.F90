! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck
! @author Slavko Brdar

#include "fckit/fctest.h"



! -----------------------------------------------------------------------------

module fcta_Field_gpu_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none

contains

subroutine module_acc_routine(view)

  implicit none
  real(4), intent(inout) :: view(:,:)

  !$acc data present(view)
  !$acc kernels
  view(1,1) = 4.
  !$acc end kernels
  !$acc end data

end subroutine module_acc_routine

end module

! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_Field_gpu,fcta_Field_gpu_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_initialize()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_host_data )
implicit none
type(atlas_Field) :: field
real(4), pointer :: view(:,:)

!!! WARNING !!! Without this interface, there is a runtime error !!!
interface
  subroutine external_acc_routine(view)
    real(4), intent(inout) :: view(:,:)
  end subroutine external_acc_routine
end interface

field = atlas_Field(kind=atlas_real(4),shape=[5,3])

call field%data(view)
view(:,:) = 0
view(1,1) = 1
call field%update_device()

!$acc data present(view)
!$acc kernels
view(1,1) = 2.
!$acc end kernels
!$acc end data

FCTEST_CHECK_EQUAL( view(1,1), 1. )
call field%update_host()
FCTEST_CHECK_EQUAL( view(1,1), 2. )

view(1,1) = 3.

call field%update_device()

write(0,*) "Calling module_acc_routine ..."
call module_acc_routine(view)
write(0,*) "Calling module_acc_routine ... done"

write(0,*) "Calling external_acc_routine ..."
call external_acc_routine(view)
write(0,*) "Calling external_acc_routine ... done"

FCTEST_CHECK_EQUAL( view(1,1), 3. )
call field%update_host()
FCTEST_CHECK_EQUAL( view(1,1), 4. )

call field%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

