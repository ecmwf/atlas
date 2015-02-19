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

call atlas_err_set_throws(.False.)
call atlas_err_set_aborts(.False.)
call atlas_err_set_backtrace(.True.)

CHECK_EQUAL(atlas_err_code(),atlas_err_cleared)

call atlas_err_success()
CHECK(atlas_noerr())
CHECK_EQUAL(atlas_err_code(),atlas_err_noerr)
CHECK_EQUAL(atlas_err_msg(),"")

call atlas_throw_exception("reason",code_location(__FILE__,__LINE__))
CHECK( atlas_err() )
CHECK_EQUAL( atlas_err_code(), atlas_err_exception )

call atlas_throw_seriousbug("reason",code_location(__FILE__,__LINE__))
CHECK_EQUAL( atlas_err_code(), atlas_err_seriousbug )

call atlas_throw_usererror("reason",code_location(__FILE__,__LINE__))
CHECK_EQUAL( atlas_err_code(), atlas_err_usererror )

call atlas_throw_outofrange("myarray",128,100)
call atlas_throw_outofrange("myarray",128,100,code_location(__FILE__,__LINE__))

CHECK_EQUAL( atlas_err_code(), atlas_err_outofrange )
! call atlas_log%warning(atlas_err_msg())
call atlas_err_clear()
CHECK_EQUAL(atlas_err_code(),atlas_err_cleared)

! call atlas_abort("I give up!!!", code_location(__FILE__,__LINE__,"test_error"))
! call atlas_abort()

! if( atlas_err() ) then
!   write(0,'(A)') atlas_err_msg()
!   call atlas_err_clear()
! endif

! call atlas_throw_exception("exception from fortran",code_location(__FILE__,__LINE__))
! call atlas_throw_notimplemented(code_location(__FILE__,__LINE__))


END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

