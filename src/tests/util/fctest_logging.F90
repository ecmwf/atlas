! (C) Copyright 2013 ECMWF.
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

module fcta_logging_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg

end module fcta_logging_fxt

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_logging,fcta_logging_fxt)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_logging )


call atlas_log%info("test_logging begin")
call atlas_log%debug("hello world")
call atlas_log%warning("World is ending")
write(msg,'(A)') "goodbye"
call atlas_log%debug(msg,newl=.False.)
call atlas_log%debug(" world")
call atlas_log%error(" oops ")
call atlas_log%info("test_logging end")
call atlas_log%panic("AAAAAGRRG")
END_TEST

! -----------------------------------------------------------------------------

!TEST( test_channel )
!type( atlas_logchannel ) :: info, debug, warning, error
!call atlas_log%channel_info%log("atlas_log%channel_info%log() called")
!info    = atlas_log%channel(ATLAS_LOG_CATEGORY_INFO)
!error   = atlas_log%channel(ATLAS_LOG_CATEGORY_ERROR)
!call info%log("info%log() called")
!call error%log("error%log() called")

!END_TEST

! -----------------------------------------------------------------------------

TEST( test_log_fortran_unit )

integer :: NULLOUT = 9
OPEN( NULLOUT, FILE='null.out' )

call atlas_log%set_fortran_unit(NULLOUT)

call atlas_log%error("error to fortran")
call atlas_log%info("info to fortran")
call atlas_log%debug("debug to fortran")
!call atlas_log%stats("stats to fortran")

CLOSE( NULLOUT )

END_TEST



! -----------------------------------------------------------------------------

END_TESTSUITE

