! (C) Copyright 1996-2016 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the logging facilities
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_logging_fixture
use atlas_module
use, intrinsic :: iso_c_binding
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

TEST( test_channel )
type( atlas_logchannel ) :: info, debug, warning, error
call atlas_log%channel_info%log("atlas_log%channel_info%log() called")
info    = atlas_log%channel(ATLAS_LOG_CATEGORY_INFO)
error   = atlas_log%channel(ATLAS_LOG_CATEGORY_ERROR)
call info%log("info%log() called")
call error%log("error%log() called")

END_TEST

! -----------------------------------------------------------------------------

TEST( test_log_fortran_unit )

integer :: NULLERR = 0
integer :: NULLOUT = 9
OPEN( NULLOUT, FILE='null.out' )

call atlas_log%channel_error%connect_fortran_unit(NULLERR)
call atlas_log%connect_fortran_unit( NULLOUT )

call atlas_log%set_prefix_stdout("[%P] info -- ")
call atlas_log%set_prefix_fortran_unit(NULLOUT, "[%5P] (%H:%M:%S) -- ")

call atlas_log%set_prefix_stdout("[%P] stdout -- ")
call atlas_log%set_prefix_stderr("[%P] stdout -- ")
call atlas_log%set_prefix_fortran_unit(NULLOUT,"[%P] NULLOUT -- ")
call atlas_log%channel_error%set_prefix_fortran_unit(NULLERR,"[%P] NULLERR -- ")

!call atlas_log%disconnect_stderr()
!call atlas_log%disconnect_stdout()
!call atlas_log%disconnect_fortran_unit(NULLOUT )

call atlas_log%error("error to fortran")
call atlas_log%info("info to fortran")
call atlas_log%debug("debug to fortran")
call atlas_log%stats("stats to fortran")

CLOSE( NULLOUT )

END_TEST



! -----------------------------------------------------------------------------

END_TESTSUITE

