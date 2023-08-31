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

TESTSUITE_WITH_FIXTURE(fctest_atlas_MultiField,fcta_MultiField_fixture)

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

    type(atlas_MultiField) :: mfield
    type(atlas_config)     :: config

    integer, parameter :: nproma  = 16;
    integer, parameter :: nlev    = 100;
    integer, parameter :: ngptot  = 2000;
    type(atlas_Config), dimension(5) :: field_configs

    config = atlas_Config()
    call config%set("type", "MultiFieldCreatorIFS");
    call config%set("ngptot", ngptot);
    call config%set("nproma", nproma);
    call config%set("nlev", nlev);
    call config%set("datatype", "real64");
    field_configs(1) = atlas_Config()
    field_configs(2) = atlas_Config()
    field_configs(3) = atlas_Config()
    field_configs(4) = atlas_Config()
    field_configs(5) = atlas_Config()
    call field_configs(1)%set("name", "temperature")
    call field_configs(2)%set("name", "pressure")
    call field_configs(3)%set("name", "density")
    call field_configs(4)%set("name", "clv")
    call field_configs(4)%set("nvar", 5)
    call field_configs(5)%set("name", "wind_u")
    call config%set("fields", field_configs)

    mfield = atlas_MultiField(config)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
