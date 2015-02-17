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

module fctest_atlas_logging_fixture
use atlas_module
use iso_c_binding
implicit none

end module fctest_atlas_logging_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_logging,fctest_atlas_logging_fixture)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_logging )
call atlas_log%info("test_logging begin")
call atlas_log%debug("hello world",lvl=0)
call atlas_log%warning("World is ending")
write(atlas_log%msg,'(A)') "goodbye"
call atlas_log%debug(endl=.False.,lvl=1)
call atlas_log%debug(" world",lvl=1)
call atlas_log%error(" oops ")
call atlas_log%info("test_logging end")
call atlas_log%panic("AAAAAGRRG")
END_TEST

! -----------------------------------------------------------------------------

TEST( test_log_fortran_unit )

integer :: NULLERR = 0
integer :: NULLOUT = 9
OPEN( NULLOUT, FILE='null.out' )

call atlas_log%connect_fortran_unit(atlas_log%cat_error, NULLERR )
call atlas_log%connect_fortran_unit(atlas_log%cat_all,   NULLOUT )

call atlas_log%error("error to fortran")
call atlas_log%info("info to fortran")
call atlas_log%debug("debug to fortran")
call atlas_log%stats("stats to fortran")

CLOSE( NULLOUT )

END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

