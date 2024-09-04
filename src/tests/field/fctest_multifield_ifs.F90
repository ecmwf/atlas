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

    type(atlas_MultiField)  :: mfield(2)  ! mfield(1): created through config
                                          ! mfield(2): created directly in atlas_MultiField constructor
    type(atlas_FieldSet)    :: fieldset_1, fieldset_2
    type(atlas_Field)       :: field
    type(atlas_config)      :: config
    integer, pointer :: fdata_int_2d(:,:)
    real(c_float),  pointer :: fdata_real32_2d(:,:)
    real(c_double), pointer :: fdata_real64_2d(:,:)
    real(c_double), pointer :: fdata_real64(:,:,:)

    integer, parameter :: nproma  = 16;
    integer, parameter :: nlev    = 0;
    integer, parameter :: ngptot  = 2000;
    integer, parameter :: nblk    = (ngptot + nproma - 1) / nproma
    type(atlas_Config), allocatable :: field_configs(:)
    integer :: i
    character(len=12), parameter, dimension(5) :: var_names = (/ &
        "temperature ", "pressure    ", "density     ", "clv         ", "wind_u      " &
    /)

    integer :: ivar
    return

    config = atlas_Config()
    call config%set("type", "MultiFieldCreatorIFS");
    call config%set("ngptot", ngptot);
    call config%set("nproma", nproma);
    call config%set("nlev", nlev);
    call config%set("datatype", "real64");
    allocate(field_configs(size(var_names)))
    do ivar = 1, size(var_names)
        field_configs(ivar) = atlas_Config()
        call field_configs(ivar)%set("name", trim(var_names(ivar)))
    end do
    !call field_configs(4)%set("nvar", 4) ! clv has four subvariables
    call config%set("fields", field_configs)
    mfield(1) = atlas_MultiField(config)

    !mfield(2) = atlas_MultiField('real64', [nproma, nlev, -1, nblk], var_names)

    do i = 1, 1
      FCTEST_CHECK_EQUAL(mfield(i)%size(), 5)

      fieldset_1 = mfield(i)%fieldset()
      FCTEST_CHECK_EQUAL(fieldset_1%size(), 5)

      fieldset_2 = atlas_FieldSet()
      call fieldset_2%add(mfield(i)%fieldset())
      field = fieldset_2%field("density")
      call field%data(fdata_real64)
      fdata_real64(1,1,1) = 2.
      call field%rename("dens")

      ! check data access directly though multifield
      call mfield(i)%data("dens", fdata_real64)
      fdata_real64(1,1,1) = 3.

      field = fieldset_1%field("dens")
      call field%data(fdata_real64)
      FCTEST_CHECK_EQUAL(fdata_real64(1,1,1), 3._c_double)
    end do

END_TEST


TEST( test_multifield_direct_constructor )
    implicit none

    type(atlas_MultiField)  :: mfield(2)
    type(atlas_FieldSet)    :: fset
    type(atlas_Field)       :: field
    type(atlas_config)      :: config
    real(c_float), pointer  :: fdata_f2d(:,:), fdata_f3d(:,:,:)
    real(c_double), pointer :: fdata_d3d(:,:,:)

    integer, parameter :: nproma  = 16;
    integer, parameter :: nlev    = 100;
    integer, parameter :: ngptot  = 2000;
    integer, parameter :: nblk    = (ngptot + nproma - 1) / nproma
    integer :: i
    character(len=12), parameter, dimension(5) :: var_names = (/ &
        "temperature ", "pressure    ", "density     ", "clv         ", "wind_u      " &
    /)
    integer :: ivar

    config = atlas_Config()

    ! surface vars
    mfield(1) = atlas_MultiField('real32', [nproma, -1, nblk], var_names, config)

    ! 3d vars
    mfield(2) = atlas_MultiField('real64', [nproma, nlev, -1, nblk], var_names, config)

    FCTEST_CHECK_EQUAL(mfield(1)%size(), 5)

    call mfield(1)%data("density", fdata_f2d)
    fdata_f2d(1,1) = 3.
    call mfield(2)%data("density", fdata_d3d)
    fdata_d3d(1,1,1) = 4.

    fset = atlas_FieldSet()
    call fset%add(mfield(1)%fieldset())
    call fset%add(mfield(2)%fieldset())

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
