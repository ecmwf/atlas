! (C) Copyright 2013 ECMWF.
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

TEST( test_fieldset_cummulative_gpu_calls_all_fields )
implicit none
type(atlas_FieldSet) :: fset
type(atlas_Field) :: field
real(4), pointer :: fview(:,:)

print *, "test_fieldset_cummulative_gpu_calls_all_fields"
fset = atlas_FieldSet()
call fset%add(atlas_Field(name="f1", kind=atlas_real(4), shape=[1,2]))
call fset%add(atlas_Field(name="f2", kind=atlas_real(4), shape=[2,2]))
call fset%add(atlas_Field(name="f3", kind=atlas_real(4), shape=[2,1]))

print *, "... by idx"
field = fset%field(2)
call field%data(fview)
fview(:,:) = 1.
fview(2,1) = 2.

call fset%set_host_needs_update(.false.)
FCTEST_CHECK_EQUAL(field%device_allocated(), .false.)
call fset%allocate_device()
FCTEST_CHECK_EQUAL(field%device_allocated(), .true.)
call fset%update_device()

!$acc data present(fview)
!$acc kernels
fview(2,1) = 5.
!$acc end kernels
!$acc end data

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%update_host()
FCTEST_CHECK_EQUAL( fview(2,1), 5. )
call fset%deallocate_device()
FCTEST_CHECK_EQUAL(field%device_allocated(), .false.)

print *, "... by name"
field = fset%field("f3")
call field%data(fview)
fview(:,:) = 3.
fview(2,1) = 4.

call fset%allocate_device()
call fset%update_device()

!$acc data present(fview)
!$acc kernels
fview(2,1) = 7.
!$acc end kernels
!$acc end data

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%update_host()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset_cummulative_gpu_calls_subset_field)
implicit none
type(atlas_FieldSet) :: fset
type(atlas_Field) :: field
real(4), pointer :: fview(:,:)
print *, "fieldset_cummulative_gpu_calls subset fields"

fset = atlas_FieldSet()
call fset%add(atlas_Field(name="f1", kind=atlas_real(4), shape=[1,2]))
call fset%add(atlas_Field(name="f2", kind=atlas_real(4), shape=[2,2]))
call fset%add(atlas_Field(name="f3", kind=atlas_real(4), shape=[2,1]))

print *, "... by idx"
field = fset%field(2)
call field%data(fview)
fview(:,:) = 1.
fview(2,1) = 2.

call fset%set_host_needs_update((/2, 3/), .false.)
call fset%allocate_device((/ 2, 3 /)) ! device allocate only for field 2 and 3
call fset%sync_host_device((/ 2, 3 /)) ! device-memcpy for field 2 and 3

!$acc data present(fview)
!$acc kernels
fview(2,1) = 5.
!$acc end kernels
!$acc end data

call fset%set_host_needs_update((/ 2 /))

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%sync_host_device()
FCTEST_CHECK_EQUAL( fview(2,1), 5. )
call fset%deallocate_device()

print *, "... by name"
field = fset%field("f3")
call field%data(fview)
fview(:,:) = 3.
fview(2,1) = 4.

call fset%allocate_device((/ "f2", "f3" /))
call fset%sync_host_device((/ "f2", "f3" /))

!$acc data present(fview)
!$acc kernels
fview(2,1) = 7.
!$acc end kernels
!$acc end data

call fset%set_host_needs_update((/ "f3" /))

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%sync_host_device()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset_cummulative_gpu_calls_not_all_field_with_sync_host_device)
implicit none
type(atlas_FieldSet) :: fset
type(atlas_Field) :: field
real(4), pointer :: fview(:,:)

print *, "test_fieldset_cummulative_gpu_calls_with_sync_host_device"
fset = atlas_FieldSet()
call fset%add(atlas_Field(name="f1", kind=atlas_real(4), shape=[1,2]))
call fset%add(atlas_Field(name="f2", kind=atlas_real(4), shape=[2,2]))
call fset%add(atlas_Field(name="f3", kind=atlas_real(4), shape=[2,1]))

print *, "... by idx"
field = fset%field(2)
call field%data(fview)
fview(:,:) = 1.
fview(2,1) = 2.

call fset%set_host_needs_update(.false.)
call fset%allocate_device((/ 2, 3 /)) ! device allocate only for field 2 and 3
call fset%set_device_needs_update((/ 2 /))
call fset%sync_host_device((/ 2 /)) ! device-memcpy for field 2 and 3

!$acc data present(fview)
!$acc kernels
fview(2,1) = 5.
!$acc end kernels
!$acc end data

call fset%set_host_needs_update((/ 2 /))

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%sync_host_device()
FCTEST_CHECK_EQUAL( fview(2,1), 5. )

call fset%deallocate_device()

print *, "... by name"
field = fset%field("f2")
call field%data(fview)
fview(:,:) = 3.
fview(2,1) = 4.

call fset%set_host_needs_update(.false.)
call fset%allocate_device((/ "f2", "f3" /)) ! device allocate only for field 2 and 3
call fset%set_device_needs_update((/ "f2" /)) ! set fields for sync_host_device
call fset%sync_host_device((/ "f2" /)) ! device-memcpy for field 'f2'

!$acc data present(fview)
!$acc kernels
fview(2,1) = 7.
!$acc end kernels
!$acc end data

call fset%set_host_needs_update((/ "f2" /))

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%sync_host_device()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

