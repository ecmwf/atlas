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

  !$acc data present(fview)
  !$acc parallel loop
  do j=1,N
    fview(1,1,j) = real(j, c_double)
    fview(2,1,j) = -3_c_double
  enddo
  !$acc end data

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
  real(c_double), pointer :: fview(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: i,j,k,l
  write(0,*) "test_field_wrapdataslice [skipped]" ! NOT DONE YET !!!
  return ! SKIP THIS TEST !!!
  allocate( existing_data(4,3,2,5) )

  existing_data = -1.

  field = atlas_Field(existing_data(:,:,1,:))
  call field%data(fview)

  !call field%allocate_device()
  !call field%update_device()

  !!$acc data present(fview)
  !!$acc parallel loop
  do i=1,4
    do j=1,3
      do k=1,2
        do l=1,5
          fview(i,j,l) = 1000*i+100*j+10*1+l
        enddo
      enddo
    enddo
  enddo
  !!$acc end data

  call field%deallocate_device()

  k=1
  do i=1,4
    do j=1,3
      do l=1,5
        FCTEST_CHECK_EQUAL(fview(i,j,l) , real(1000*i+100*j+10*k+l,c_double) )
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

  !$acc data present(fview)
  !$acc parallel loop
  do i=1,Ni
    do j=1,Nj
      do k=1,Nk
        fview(i,j,k) = (mod(k,3) == 0 )
      enddo
    enddo
  enddo
  !$acc end data

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

