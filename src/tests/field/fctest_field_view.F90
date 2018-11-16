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

module fcta_Field_view_fixture
use atlas_module
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env
implicit none
end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Field_view,fcta_Field_view_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_field_data_view)
use omp_lib

implicit none

integer, parameter :: lon = 20, lev = 16, blocks = 138
type(atlas_Field) :: field
real(kind=REAL64), pointer :: global(:,:,:), local(:,:)
integer :: i, j, k, pid, num_threads

field = atlas_Field(kind=atlas_real(kind=REAL64), shape=[lon, lev, blocks])
FCTEST_CHECK_EQUAL( field%rank()   , 3  )
FCTEST_CHECK_EQUAL( field%shape(1) , lon  )
FCTEST_CHECK_EQUAL( field%shape(2) , lev  )
FCTEST_CHECK_EQUAL( field%shape(3) , blocks  )

!$OMP PARALLEL PRIVATE(k, global, local) SHARED(field)
num_threads = omp_get_num_threads()

!$OMP DO SCHEDULE(DYNAMIC,1)
do k=1, blocks
  pid = omp_get_thread_num()

  call field%data(global)
  local => global(:,:,k)

  local(:,:) = real(pid, REAL64)
end do
!$OMP END DO
!$OMP END PARALLEL

call field%data(global)

do k=1, blocks
  do j=1, lev
    do i=1, lon
      FCTEST_CHECK( global(i,j,k) >= 0.0 )
      FCTEST_CHECK( global(i,j,k) <= REAL(num_threads, REAL64) )
    end do
  end do
end do

END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

