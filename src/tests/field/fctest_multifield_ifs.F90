! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck
! @author Slavko Brdar

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_MultiField_fixture
use atlas_module
use atlas_multifield_module
use, intrinsic :: iso_c_binding
implicit none
end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_MultiField, fcta_MultiField_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_multifield )
    implicit none

    type(atlas_MultiField)  :: multifield
    type(atlas_FieldSet)    :: fieldset_1, fieldset_2
    type(atlas_Field)       :: field
    type(atlas_config)      :: config
    real(c_double), pointer :: view(:,:,:)

    integer, parameter :: nvar = 5;
    integer, parameter :: nproma  = 16;
    integer, parameter :: nlev    = 100;
    integer, parameter :: ngptot  = 2000;
    type(atlas_Config), dimension(5) :: field_configs
    character(len=12), parameter, dimension(nvar) :: var_names = (/ &
        "temperature ", "pressure    ", "density     ", "clv         ", "wind_u      " &
    /)

    integer :: ivar

    config = atlas_Config()
    call config%set("type", "MultiFieldCreatorIFS");
    call config%set("ngptot", ngptot);
    call config%set("nproma", nproma);
    call config%set("nlev", nlev);
    call config%set("datatype", "real64");
    do ivar = 1, 5
        field_configs(ivar) = atlas_Config()
        call field_configs(ivar)%set("name", trim(var_names(ivar)))
    end do
    call field_configs(4)%set("nvar", 5) ! clv has five subvariables
    call config%set("fields", field_configs)

    multifield = atlas_MultiField(config)
    FCTEST_CHECK_EQUAL(multifield%size(), 5)

    fieldset_1 = multifield%fieldset()
    FCTEST_CHECK_EQUAL(fieldset_1%size(), 5)

    fieldset_2 = atlas_FieldSet()
    call fieldset_2%add(multifield%fieldset())
    field = fieldset_2%field("density")
    call field%data(view)
    view(1,1,1) = 2

    field = fieldset_1%field("density")
    call field%data(view)
    FCTEST_CHECK_EQUAL(view(1,1,1), 2)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
