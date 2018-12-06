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
#include "atlas/atlas_f.h"

! -----------------------------------------------------------------------------

module fcta_Metadata_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_Metadata_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Metadata,fcta_Metadata_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use fckit_main_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_metadata )
use fckit_log_module
use fckit_c_interop_module
implicit none

  integer(c_int) :: intval
  logical :: true, false
  real(c_float) :: real32
  real(c_double) :: real64
  character(len=:), allocatable :: string
  integer(c_int), allocatable :: arr_int32(:)
  real(c_float), allocatable :: arr_real32(:)
  type(atlas_Metadata) :: metadata
  type(fckit_logchannel) :: info

  write(*,*) "test_metadata starting"

  metadata = atlas_Metadata()

  write(0,*) "metadata%c_ptr() = ", c_ptr_to_loc(metadata%CPTR_PGIBUG_A)

  call metadata%set("true",.True.)
  call metadata%set("false",.False.)
  call metadata%set("int",20)
  call metadata%set("real32", real(0.1,kind=c_float)   )
  call metadata%set("real64", real(0.2,kind=c_double) )
  call metadata%set("string", "hello world")
  call metadata%set("arr_int32", (/1,2,3/))
  call metadata%set("arr_int64", (/1_c_long,2_c_long,3_c_long/))
  call metadata%set("arr_real32", (/1.1_c_float,2.1_c_float,3.7_c_float/))
  call metadata%set("arr_real64", (/1.1_c_double,2.1_c_double,3.7_c_double/))

  call metadata%get("true",true)
  call metadata%get("false",false)
  call metadata%get("int",intval)
  call metadata%get("real32",real32)
  call metadata%get("real64",real64)
  call metadata%get("string",string)
  call metadata%get("arr_int64",arr_int32)
  call metadata%get("arr_real64",arr_real32)

  info = fckit_log%info_channel()
  call metadata%print(info)
  call fckit_log%info("",newl=.true.)
  call fckit_log%info(metadata%json())

  CHECK( true  .eqv. .True.  )
  CHECK( false .eqv. .False. )

  FCTEST_CHECK_EQUAL( intval, 20 )
  FCTEST_CHECK_CLOSE( real32, real(0.1,kind=c_float), real(0.,kind=c_float) )
  FCTEST_CHECK_CLOSE( real64, real(0.2,kind=c_double), real(0.,kind=c_double) )
  FCTEST_CHECK_EQUAL( string, "hello world" )
  FCTEST_CHECK_EQUAL( arr_int32, (/1,2,3/) )
  FCTEST_CHECK_EQUAL( arr_real32, (/1.1_c_float,2.1_c_float,3.7_c_float/) )

  call metadata%final()

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

