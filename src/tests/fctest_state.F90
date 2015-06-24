! (C) Copyright 1996-2014 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_State_Fixture
use atlas_module
use iso_c_binding
implicit none

contains

end module fctest_atlas_State_Fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_State,fctest_atlas_State_Fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_state_fields )
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
write(atlas_log%msg,'(A,I0,A)') "The state contains ",state%nb_fields()," fields."; call atlas_log%info()

! Check if wind field exists
if( .not. state%has_field("wind") ) then
  write(atlas_log%msg,'(A)') "The state does not contain the wind field"; call atlas_log%info()
endif

! Print existing fields info
write(atlas_log%msg,'(A)') "The state contains the fields:"; call atlas_log%info()
do jfield=1,state%nb_fields()
  field = state%field(jfield)
  write(atlas_log%msg,'(2A)')   "  - ",field%name();                   call atlas_log%info()
  write(atlas_log%msg,'(2A)')   "        kind = ",field%data_type();   call atlas_log%info()
  write(atlas_log%msg,'(A,I0)') "        size = ",field%size();        call atlas_log%info()
  write(atlas_log%msg,*)         "       shape =",field%shape();       call atlas_log%info()
enddo

! Get fields out of the state
temperature_field = state%field("temperature")
pressure_field = state%field("pressure")

! If you want to edit the field, access the data
call temperature_field%access_data(temperature)
call pressure_field%access_data(pressure)

temperature(:,:) = 273.15_c_double
pressure(:,:,:) = 10000._c_double

! Set some field metadata
metadata = temperature_field%metadata()
call metadata%set("unit","Kelvin")
call metadata%set("iteration",1)
call metadata%set("grib_param_id","T")

call state%remove_field("pressure")
FCTEST_CHECK(state%has_field("pressure") .eqv. .False.)

! Delete the state
call atlas_delete(state)

END_TEST

! -----------------------------------------------------------------------------

TEST( test_state_grids )
type(atlas_State) :: state
type(atlas_Field) :: field
type(atlas_ReducedGrid) :: grid

! Create a new state
state = atlas_State()

! Create a new grid inside the state
call state%add( atlas_ReducedGrid("oct.N48") )

! Check how many fields we have
write(atlas_log%msg,'(A,I0,A)') "The state contains ",state%nb_grids()," grid."; call atlas_log%info()

grid = state%grid()

write(atlas_log%msg,*) "The grid has ",grid%npts()," points"; call atlas_log%info()

! Delete the state
call atlas_delete(state)

END_TEST


! -----------------------------------------------------------------------------

TEST( test_state_meshes )
type(atlas_State) :: state
type(atlas_Field) :: field
type(atlas_mesh) :: mesh
type(atlas_FunctionSpace) :: nodes

! Create a new state
state = atlas_State()

! Create a new grid and mesh inside the state
call state%add( atlas_ReducedGrid("oct.N48") )
call state%add( atlas_generate_mesh( state%grid() ) )

! Check how many fields we have
write(atlas_log%msg,'(A,I0,A)') "The state contains ",state%nb_meshes()," mesh."; call atlas_log%info()

mesh = state%mesh()
nodes = mesh%function_space("nodes")

write(atlas_log%msg,*) "The mesh has ",nodes%dof()," points"; call atlas_log%info()

! Delete the state
call atlas_delete(state)

END_TEST

! -----------------------------------------------------------------------------


TEST( test_state_factory )
type(atlas_State) :: state
type(atlas_Field) :: field
type(atlas_mesh) :: mesh
type(atlas_FunctionSpace) :: nodes

! Create a new state
state = atlas_State()

! Delete the state
call atlas_delete(state)

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

