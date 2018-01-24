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

! -----------------------------------------------------------------------------

module fctest_atlas_State_Fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none

  character(len=1024) :: msg

contains

end module fctest_atlas_State_Fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_State,fctest_atlas_State_Fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_state_fields )
#if 1
type(atlas_State) :: state
type(atlas_Field) :: temperature_field
type(atlas_Field) :: pressure_field
type(atlas_Field) :: field
type(atlas_Metadata) :: metadata

integer :: jfield

real(c_double), pointer :: temperature(:,:)
real(c_float), pointer :: pressure(:,:,:)

! Create a new state
state = atlas_State()

! Create a new fields inside the state
call state%add( atlas_Field(shape=[20,10],   kind=atlas_real(c_double), name="temperature") )
call state%add( atlas_Field(shape=[20,10,3], kind=atlas_real(c_float),  name="pressure"   ) )

! Check how many fields we have
write(msg,'(A,I0,A)') "The state contains ",state%size()," fields."; call atlas_log%info(msg)
FCTEST_CHECK_EQUAL( state%size(), 2 )

! Check if wind field exists
if( .not. state%has("wind") ) then
  write(msg,'(A)') "The state does not contain the wind field"; call atlas_log%info(msg)
endif
FCTEST_CHECK( .not. state%has("wind") )

! Print existing fields info
write(msg,'(A)') "The state contains the fields:"; call atlas_log%info(msg)
do jfield=1,state%size()
  field = state%field(jfield)
  write(msg,'(2A)')   "  - ",field%name();                   call atlas_log%info(msg)
  write(msg,'(2A)')   "        kind = ",field%datatype();    call atlas_log%info(msg)
  write(msg,'(A,I0)') "        size = ",field%size();        call atlas_log%info(msg)
  write(msg,*)         "       shape =",field%shape();       call atlas_log%info(msg)
enddo

! Get fields out of the state
temperature_field = state%field("temperature")
pressure_field = state%field("pressure")

! If you want to edit the field, access the data
call temperature_field%data(temperature)
call pressure_field%data(pressure)

temperature(:,:) = 273.15_c_double
pressure(:,:,:) = 10000._c_double

! Set some field metadata
metadata = temperature_field%metadata()
call metadata%set("unit","Kelvin")
call metadata%set("iteration",1)
call metadata%set("grib_param_id","T")

call state%remove("pressure")
FCTEST_CHECK(state%has("pressure") .eqv. .False.)

! Delete the state
call state%final()
#endif
END_TEST

! -----------------------------------------------------------------------------


TEST( test_state_factory )
#if 1
type(atlas_State) :: state

! Create a new state
state = atlas_State()

! Delete the state
call state%final()
#endif
END_TEST

! -----------------------------------------------------------------------------

TEST( test_state_metadata )
#if 1
type(atlas_State) :: state
type(atlas_Metadata) :: state_metadata

! Create a new state
state = atlas_State()

! Access metadata
state_metadata = state%metadata()

call state_metadata%set("integer",1)

! Delete the state
call state%final()
#endif
END_TEST

! -----------------------------------------------------------------------------
END_TESTSUITE

