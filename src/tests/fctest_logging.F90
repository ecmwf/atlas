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
call logger%info("test_logging begin")
call logger%debug("hello world",lvl=0)
call logger%warning("World is ending")
write(logger%message,'(A)') "goodbye"
call logger%debug(endl=.False.,lvl=1)
call logger%debug(" world",lvl=1)
call logger%error(" oops ")
call logger%info("test_logging end")
call logger%panic("AAAAAGRRG")
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

