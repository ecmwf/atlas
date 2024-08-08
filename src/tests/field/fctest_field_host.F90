! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"
#include "atlas/atlas_f.h"

! -----------------------------------------------------------------------------

module fcta_Field_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none

contains

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_Field,fcta_Field_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_initialize()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_host_data )
type(atlas_Field) :: field
real(8), pointer :: view(:,:)

field = atlas_Field(kind=atlas_real(8),shape=[10,5])

call field%data(view)

FCTEST_CHECK( .not. field%host_needs_update() )
#if ! ATLAS_HAVE_CUDA
FCTEST_CHECK( .not. field%device_needs_update() )
#endif

call field%update_device()
FCTEST_CHECK( .not. field%device_needs_update() )

call field%final()
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

