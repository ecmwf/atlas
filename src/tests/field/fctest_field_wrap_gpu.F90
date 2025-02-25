! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_Field_wrap_device_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Field_wrap_device,fcta_Field_wrap_device_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

subroutine kernel_1(view, N)
  real(c_double) :: view(:,:,:)
  integer, intent(in) :: N
  !$acc data present(view)
  !$acc parallel loop
  do j=1,N
    view(1,1,j) = real(j, c_double)
    view(2,1,j) = -3_c_double
  enddo
  !$acc end data
end subroutine kernel_1


subroutine kernel_2(view, Ni, Nj, Nl)
  ! The OpenACC standard states that a deviceptr (here 'view')
  ! shall not have the POINTER attribute, hence this subroutine
  ! is required rather than have the code inlined at call site.
  real(c_double) :: view(:,:,:)
  integer, intent(in) :: Ni, Nj, Nl
  !$acc data deviceptr(view)
  !$acc parallel loop
  do i = 1, Ni
    do j = 1, Nj
      do l = 1, Nl
        view(i,j,l) = -view(i,j,l)
      enddo
    enddo
  enddo
  !$acc end data
end subroutine kernel_2


subroutine kernel_3(view, Ni, Nj, Nk)
  logical :: view(:,:,:)
  integer, intent(in) :: Ni, Nj, Nk
  !$acc data present(view)
  !$acc parallel loop
  do i=1,Ni
    do j=1,Nj
      do k=1,Nk
        view(i,j,k) = (mod(k,3) == 0 )
      enddo
    enddo
  enddo
  !$acc end data
end subroutine kernel_3


TEST( test_field_wrapdata )
implicit none

  real(c_double), allocatable :: existing_data(:,:,:)
  real(c_double), pointer :: fview(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: N=20
  integer(c_int) :: j
  write(0,*) "test_field_wrapdata"
  allocate( existing_data(2,10,N) )

  existing_data(:,:,:) = -1._c_double

  field = atlas_Field("wrapped", existing_data)
  call field%data(fview)

  do j=1,N
    fview(1,1,j) = -2._c_double
  enddo

  call field%allocate_device()
  call field%update_device()

  call kernel_1(fview, N)

  j = N/2
  FCTEST_CHECK_EQUAL( existing_data(1,1,j), -2._c_double )
  FCTEST_CHECK_EQUAL( existing_data(2,1,j), -1._c_double)

  call field%update_host()
  call field%deallocate_device()

  call field%final()

  ! Existing data is not deleted after field%final()
  FCTEST_CHECK_EQUAL( existing_data(1,1,j), real(j, c_double) )
  FCTEST_CHECK_EQUAL( existing_data(2,1,j), -3._c_double )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_wrapdataslice )
implicit none

  real(c_double), allocatable :: existing_data(:,:,:,:)
  real(c_double), pointer :: dview(:,:,:), hview(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: i, j, k, l
  integer(c_int), parameter :: Ni = 8, Nj = 2, Nk = 3, Nl = 4, k_idx = 2
  real(c_double) :: a
  write(0,*) "test_field_wrapdataslice"
  allocate( existing_data(Ni, Nj, Nk, Nl) )

  do i = 1, Ni
    do j = 1, Nj
      do k = 1, Nk
        do l = 1, Nl
          existing_data(i,j,k,l) = real(1000 * i+100 * j + 10*k + l, c_double)
        enddo
      enddo
    enddo
  enddo

  field = atlas_Field(existing_data(:,:,k_idx,:))
  FCTEST_CHECK_EQUAL(field%shape() , ([Ni, Nj, Nl]))
  FCTEST_CHECK_EQUAL(field%strides() , ([1, Ni, Ni*Nj*Nk]))

  call field%allocate_device()
  call field%update_device()

  call field%data(hview)
  call field%device_data(dview)

  FCTEST_CHECK_EQUAL(shape(hview), ([Ni, Nj, Nl]))
  FCTEST_CHECK_EQUAL(shape(dview), ([Ni, Nj, Nl]))

  call kernel_2(dview, Ni, Nj, Nl)

  call field%update_host()
  call field%deallocate_device()

  do i = 1, Ni
    do j = 1, Nj
      do k = 1, Nk
        a = 1._c_double;
        if (k == k_idx) a = -1._c_double;
        do l = 1, Nl
          FCTEST_CHECK_EQUAL(existing_data(i, j, k, l), a * real(1000*i + 100*j + 10*k + l, c_double))
        enddo
      enddo
    enddo
  enddo

  call field%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_wrap_logical)
implicit none

  logical, allocatable :: existing_data(:,:,:)
  logical, pointer :: fview(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: Ni=1
  integer(c_int) :: Nj=1
  integer(c_int) :: Nk=6
  integer(c_int) :: i,j,k
  write(0,*) "test_field_wrap_logical"
  allocate( existing_data(Ni,Nj,Nk) )

  do i=1,Ni
    do j=1,Nj
      do k=1,Nk
        existing_data(i,j,k) = (mod(k,2) == 0 )
      enddo
    enddo
  enddo

  ! Work with fields from here
  field = atlas_Field("wrapped",existing_data)
  call field%data(fview)

  call field%allocate_device()
  call field%update_device()

  call kernel_3(fview, Ni, Nj, Nk)

  FCTEST_CHECK_EQUAL( fview(1,1,1), .false. )
  FCTEST_CHECK_EQUAL( fview(1,1,2), .true. )

  call field%update_host()
  call field%deallocate_device()

  call field%final()
  ! ... until here

  ! Existing data is not deleted
  do i=1,Ni
    do j=1,Nj
      do k=1,Nk
        FCTEST_CHECK_EQUAL( fview(i,j,k), (mod(k,3) == 0) )
      enddo
    enddo
  enddo
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

