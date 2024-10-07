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
  use, intrinsic :: iso_c_binding, only: c_ptr
  implicit none

  !!! WARNING !!! Without this interface, there is a runtime error !!!
  interface
    subroutine external_acc_routine(view)
      real(4), pointer, intent(inout) :: view(:,:)
    end subroutine external_acc_routine
    subroutine external_acc_routine_devptr(dview)
      real(4), pointer, intent(inout) :: dview(:,:)
    end subroutine external_acc_routine_devptr
  end interface

contains

  subroutine module_acc_routine(view)
    implicit none
    real(4), intent(inout) :: view(:,:)

    !$acc kernels present(view)
    view(1,1) = 4.
    !$acc end kernels

  end subroutine module_acc_routine

! -----------------------------------------------------------------------------

  subroutine check_field(field, memory_mapped)
    use fctest, only: fce
    implicit none

    type(atlas_Field), intent(inout) :: field
    real(4), pointer :: view(:,:)
    real(4), pointer :: dview(:,:)
    logical, intent(in) :: memory_mapped

    if (field%name() == '') call field%rename("field_no-pinning")

    call field%data(view)
    view(:,:) = 0
    view(1,1) = 1
    call field%allocate_device()
    call field%update_device()

    !$acc kernels present(view)
        view(1,1) = 2.
    !$acc end kernels

    if (.not. memory_mapped) FCTEST_CHECK_EQUAL( view(1,1), 1. )
    call field%update_host()
    FCTEST_CHECK_EQUAL( view(1,1), 2. )

    call field%device_data(dview)
    !$acc kernels deviceptr(dview)
        view(1,1) = 3.
    !$acc end kernels

    if (.not. memory_mapped) FCTEST_CHECK_EQUAL( view(1,1), 2. )
    call field%update_host()
    FCTEST_CHECK_EQUAL( view(1,1), 3. )

    print *, "Check module_acc_routine on ", field%name()
    call module_acc_routine(view)

    if (.not. memory_mapped) FCTEST_CHECK_EQUAL( view(1,1), 3. )
    call field%update_host()
    FCTEST_CHECK_EQUAL( view(1,1), 4. )

    print *, "Check external_acc_routine on ", field%name()
    call external_acc_routine(view)

    if (.not. memory_mapped) FCTEST_CHECK_EQUAL( view(1,1), 4. )
    call field%update_host()
    FCTEST_CHECK_EQUAL( view(1,1), 5. )

    call field%deallocate_device()
  end subroutine

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

TEST( test_host_data_with_memory_pinned_mapped )
implicit none
type(atlas_Field) :: field
real(4), pointer :: view(:,:)
type(atlas_Config) :: options

! host memory pinning with mapped device memory
options = atlas_Config()
call options%set("host_memory_pinned", .true.)
call options%set("host_memory_mapped", .true.)
field = atlas_Field(name="field_pinned-mapped", kind=atlas_real(4), shape=[5,3], options=options)
call check_field(field, memory_mapped = .true.)
call field%final()

! host memory pinning, no field name
call options%set("host_memory_pinned", .true.)
call options%set("host_memory_mapped", .false.)
field = atlas_Field(kind=atlas_real(4), shape=[5,3], options=options)
call check_field(field, memory_mapped = .false.)
call field%final()

! memory no pinning
call options%set("host_memory_pinned", .false.)
field = atlas_Field(kind=atlas_real(4), shape=[5,3], options=options)
call check_field(field, memory_mapped = .false.)
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
if (ATLAS_HAVE_GRIDTOOLS_STORAGE == 0) then
  FCTEST_CHECK_EQUAL(field%device_allocated(), .false.)
endif
call fset%allocate_device()
FCTEST_CHECK_EQUAL(field%device_allocated(), .true.)
call fset%update_device()

!$acc kernels present(fview)
fview(2,1) = 5.
!$acc end kernels

FCTEST_CHECK_EQUAL( fview(2,1), 2. )
call fset%update_host()
FCTEST_CHECK_EQUAL( fview(2,1), 5. )
call fset%deallocate_device()
if (ATLAS_HAVE_GRIDTOOLS_STORAGE == 0) then
  FCTEST_CHECK_EQUAL(field%device_allocated(), .false.)
endif
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

END_TESTSUITE

