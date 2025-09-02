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

module fcta_FieldSet_gpu_fxt
  use atlas_module
  use, intrinsic :: iso_c_binding, only: c_ptr
  implicit none

contains

  subroutine kernel_host_ptr(view)
    implicit none
    real(4), intent(inout) :: view(:,:)
    !$acc data present(view)
    !$acc kernels present(view)
    view(2,1) = 4.
    !$acc end kernels
    !$acc end data
  end subroutine kernel_host_ptr

  subroutine kernel_device_ptr(dview)
    implicit none
    real(4) :: dview(:,:)
    !$acc data deviceptr(dview)
    !$acc kernels
    dview(2,1) = 5.
    !$acc end kernels
    !$acc end data
  end subroutine kernel_device_ptr

! -----------------------------------------------------------------------------

end module

! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_FieldSet_gpu,fcta_FieldSet_gpu_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_initialize()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_fieldset_cummulative_gpu_calls_all_fields )
implicit none
type(atlas_FieldSet) :: fset
type(atlas_Field) :: field
real(4), pointer :: fview(:,:)
real(4), pointer :: fview_dev(:,:)

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
call fset%allocate_device()
FCTEST_CHECK_EQUAL(field%device_allocated(), .true.)
call fset%update_device()

!$acc kernels present(fview)
fview(2,1) = 3.
!$acc end kernels

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%update_host()
FCTEST_CHECK_EQUAL( fview(2,1), 3. )

call kernel_host_ptr(fview)
FCTEST_CHECK_EQUAL( fview(2,1), 3. )
call fset%update_host()
FCTEST_CHECK_EQUAL( fview(2,1), 4. )

call field%device_data(fview_dev)
call kernel_device_ptr(fview_dev)
FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%update_host()
FCTEST_CHECK_EQUAL( fview(2,1), 5. )

call fset%deallocate_device()

print *, "... by name"
field = fset%field("f3")
call field%data(fview)
fview(:,:) = 3.
fview(2,1) = 4.

call fset%allocate_device()
call fset%update_device()

!$acc kernels present(fview)
fview(2,1) = 7.
!$acc end kernels

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

!$acc kernels present(fview)
fview(2,1) = 5.
!$acc end kernels

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

!$acc kernels present(fview)
fview(2,1) = 7.
!$acc end kernels

call fset%set_host_needs_update((/ "f3" /))

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%sync_host_device()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset_cummulative_gpu_calls_with_sync_host_device)
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

!$acc kernels present(fview)
fview(2,1) = 5.
!$acc end kernels

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

!$acc kernels present(fview)
fview(2,1) = 7.
!$acc end kernels

call fset%set_host_needs_update((/ "f2" /))

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%sync_host_device()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset_cummulative_gpu_calls_subset_field_with_sync_host_sync_device)
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
call fset%sync_device((/ 2, 3 /)) ! device-memcpy for field 2 and 3

!$acc kernels present(fview)
fview(2,1) = 5.
!$acc end kernels

call fset%set_host_needs_update((/ 2 /))

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%sync_host()
FCTEST_CHECK_EQUAL( fview(2,1), 5. )
call fset%deallocate_device()

print *, "... by name"
field = fset%field("f3")
call field%data(fview)
fview(:,:) = 3.
fview(2,1) = 4.

call fset%allocate_device((/ "f2", "f3" /))
call fset%sync_device((/ "f2", "f3" /))

!$acc kernels present(fview)
fview(2,1) = 7.
!$acc end kernels

call fset%set_host_needs_update((/ "f3" /))

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%sync_host()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset_cummulative_gpu_calls_with_sync_host_sync_device)
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
call fset%sync_device((/ 2 /)) ! device-memcpy for field 2 and 3

!$acc kernels present(fview)
fview(2,1) = 5.
!$acc end kernels

call fset%set_host_needs_update((/ 2 /))

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%sync_host()
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
call fset%sync_device((/ "f2" /)) ! device-memcpy for field 'f2'

!$acc kernels present(fview)
fview(2,1) = 7.
!$acc end kernels

call fset%set_host_needs_update((/ "f2" /))

FCTEST_CHECK_EQUAL( fview(2,1), 4. )
call fset%sync_host()
FCTEST_CHECK_EQUAL( fview(2,1), 7. )

call fset%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

