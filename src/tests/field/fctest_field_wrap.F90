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

module fcta_Field_wrap_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Field_wrap,fcta_Field_wrap_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_field_wrapdata)
implicit none

  real(c_double), allocatable :: existing_data(:,:,:)
  real(c_double), pointer :: data(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: N=20
  integer(c_int) :: j
  write(0,*) "test_field_wrapdata"
  allocate( existing_data(2,10,N) )
  write(0,*) "line ", __LINE__

  ! Work with fields from here
  field = atlas_Field("wrapped",existing_data)

  FCTEST_CHECK_EQUAL( field%rank()   , 3  )

  FCTEST_CHECK_EQUAL( field%size()   , 2*10*N )
  FCTEST_CHECK_EQUAL( field%shape(1) , 2  )
  FCTEST_CHECK_EQUAL( field%shape(2) , 10 )
  FCTEST_CHECK_EQUAL( field%shape(3) , N  )

  FCTEST_CHECK_EQUAL( field%owners() , 1  )

  call field%data(data)

  do j=1,N
    data(1,1,j) = j
  enddo

  call field%final()
  ! ... until here

  ! Existing data is not deleted
  do j=1,N
    FCTEST_CHECK_EQUAL( existing_data(1,1,j), real(j,c_double) )
  enddo
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_wrapdataslice)
implicit none

  real(c_double), allocatable :: existing_data(:,:,:,:)
  real(c_double), pointer :: data(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: i,j,k,l
  write(0,*) "test_field_wrapdataslice"
  allocate( existing_data(4,3,2,5) )

  do i=1,4
    do j=1,3
      do k=1,2
        do l=1,5
          existing_data(i,j,k,l) = 1000*i+100*j+10*k+l
        enddo
      enddo
    enddo
  enddo

  field = atlas_Field(existing_data(:,:,1,:))
  FCTEST_CHECK_EQUAL( field%rank()   , 3  )
  FCTEST_CHECK_EQUAL( field%size()   , 4*3*5 )
  FCTEST_CHECK_EQUAL( field%shape(1) , 4  )
  FCTEST_CHECK_EQUAL( field%shape(2) , 3  )
  FCTEST_CHECK_EQUAL( field%shape(3) , 5  )

  call field%data(data)

  write(0,*) "Shape of field = ",shape(data)

  k=1
  do i=1,4
    do j=1,3
      do l=1,5
        FCTEST_CHECK_EQUAL(data(i,j,l) , real(1000*i+100*j+10*k+l,c_double) )
      enddo
    enddo
  enddo

  call field%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_wrap_logical)
implicit none

  logical(c_bool), allocatable :: existing_data(:,:,:)
  logical(c_bool), pointer :: data(:,:,:)
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

  FCTEST_CHECK_EQUAL( field%rank()   , 3  )

  FCTEST_CHECK_EQUAL( field%size()   , Ni*Nj*Nk )
  FCTEST_CHECK_EQUAL( field%shape(1) , Ni  )
  FCTEST_CHECK_EQUAL( field%shape(2) , Nj )
  FCTEST_CHECK_EQUAL( field%shape(3) , Nk  )

  FCTEST_CHECK_EQUAL( field%owners() , 1  )

  call field%data(data)

  do i=1,Ni
    do j=1,Nj
      do k=1,Nk
        FCTEST_CHECK_EQUAL( logical(data(i,j,k)), (mod(k,2) == 0) )
        data(i,j,k) = (mod(k,3) == 0 )
      enddo
    enddo
  enddo

  call field%final()
  ! ... until here

  ! Existing data is not deleted
  do i=1,Ni
    do j=1,Nj
      do k=1,Nk
        FCTEST_CHECK_EQUAL( logical(data(i,j,k)), (mod(k,3) == 0) )
      enddo
    enddo
  enddo
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

