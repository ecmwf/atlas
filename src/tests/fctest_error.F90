! (C) Copyright 1996-2014 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_error_fixture
use atlas_module
use iso_c_binding
implicit none

end module fctest_atlas_error_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_error,fctest_atlas_error_fixture)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_error )
integer :: err_code
character(len=:), allocatable :: err_msg

err_code = atlas_err_code()
err_msg  = atlas_err_msg()

write(0,*) "err_code = ",err_code
write(0,*) "err_msg = ",err_msg

! call atlas_abort("I give up!!!", code_location(__FILE__,__LINE__,"test_error"))
! call atlas_abort()
call atlas_throw_exception("exception from fortran")
! call atlas_throw_exception("exception from fortran",code_location(__FILE__,__LINE__))
! call atlas_throw_notimplemented(code_location(__FILE__,__LINE__))

if( atlas_err() ) then 
  write(0,'(A)') atlas_err_msg()
  call atlas_err_clear()
endif
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

