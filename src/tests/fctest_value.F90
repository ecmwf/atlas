! (C) Copyright 1996-2015 ECMWF.
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

module fctest_atlas_value_fixture
use atlas_module
implicit none

end module fctest_atlas_value_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_value,fctest_atlas_value_fixture)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_value )
  integer(c_int) :: intval
  integer(c_long) :: longval
  real(c_float) :: floatval
  real(c_double) :: doubleval
  character(len=:), allocatable :: stringval

  type(atlas_Value) :: value

  value = atlas_Value(10_c_int)
  call value%get(intval)
  FCTEST_CHECK_EQUAL( intval , 10_c_int)
  call atlas_delete(value)

  value = atlas_Value(5000000000_c_long)
  call value%get(longval)
  FCTEST_CHECK_EQUAL(longval, 5000000000_c_long)
  call atlas_delete(value)

  value = atlas_Value(0.123456_c_float)
  call value%get(floatval)
  FCTEST_CHECK_EQUAL(floatval, 0.123456_c_float )
  call atlas_delete(value)

  value = atlas_Value(0.1234567890123456789_c_double)
  call value%get(doubleval)
  FCTEST_CHECK_EQUAL(doubleval, 0.1234567890123456789_c_double )
  call atlas_delete(value)

  value = atlas_Value("hello world")
  call value%get(stringval)
  FCTEST_CHECK_EQUAL(stringval, "hello world")
  call atlas_delete(value)

END_TEST

TEST( test_value_of_array )
  integer(c_int), allocatable :: intval(:)
  integer(c_long), allocatable :: longval(:)
  real(c_float), allocatable :: floatval(:)
  real(c_double), allocatable :: doubleval(:)

  type(atlas_Value) :: value

  value = atlas_Value((/10_c_int,11_c_int,12_c_int/))
  call value%get(intval)
  FCTEST_CHECK_EQUAL( intval , (/10_c_int,11_c_int,12_c_int/))
  call atlas_delete(value)

  value = atlas_Value((/5000000000_c_long,5000000001_c_long,5000000002_c_long/))
  call value%get(longval)
  FCTEST_CHECK_EQUAL(longval, (/5000000000_c_long,5000000001_c_long,5000000002_c_long/))
  call atlas_delete(value)

  value = atlas_Value((/0.123456_c_float,1.123456_c_float,2.123456_c_float/))
  call value%get(floatval)
  FCTEST_CHECK_EQUAL(floatval, (/0.123456_c_float,1.123456_c_float,2.123456_c_float/) )
  call atlas_delete(value)

  value = atlas_Value((/0.1234567890123456789_c_double,0.1234567890123456789_c_double,0.1234567890123456789_c_double/))
  call value%get(doubleval)
  FCTEST_CHECK_EQUAL(doubleval, (/0.1234567890123456789_c_double,0.1234567890123456789_c_double,0.1234567890123456789_c_double/) )
  call atlas_delete(value)

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

