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

module fctest_atlas_resource_fixture
use atlas_module
implicit none

end module fctest_atlas_resource_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_resource,fctest_atlas_resource_fixture)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_resource )
  integer(c_int) :: intval
  integer(c_long) :: longval
  real(c_float) :: floatval
  real(c_double) :: doubleval
  character(len=30) :: stringval

  call resource("-integer",0_c_int,intval)
  FCTEST_CHECK_EQUAL(intval, 10_c_int)
  write(0,*) "integer = ",intval

  call resource("-long",0_c_long,longval)
  write(0,*) "long = ",longval
  FCTEST_CHECK_EQUAL(longval, 5000000000_c_long)

  call resource("-float",0._c_float,floatval)
  FCTEST_CHECK_EQUAL(floatval, 0.123456_c_float )
  write(0,*) "float = ",floatval

  call resource("-double",0._c_double,doubleval)
  FCTEST_CHECK_EQUAL(doubleval, 0.1234567890123456789_c_double )
  write(0,*) "double = ",doubleval

  call resource("-string","",stringval)
  FCTEST_CHECK_EQUAL(stringval, "hello world")
  write(0,*) "string = ",stringval

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

