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
use pluto_module, only : pluto
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

    type(atlas_MultiField)  :: mfield_1, mfield_2
    type(atlas_FieldSet)    :: fieldset_1, fieldset_2
    type(atlas_Field)       :: field
    type(atlas_config)      :: config
    integer, pointer :: fdata_int_2d(:,:)
    real(c_float), pointer  :: fdata_f2d(:,:), fdata_f3d(:,:,:)
    real(c_double), pointer :: fdata_d3d(:,:,:)

    integer, parameter :: nproma  = 16;
    integer, parameter :: nlev    = 1;
    integer, parameter :: ngptot  = 2000;
    integer, parameter :: nblk    = (ngptot + nproma - 1) / nproma
    type(atlas_Config), allocatable :: field_configs(:)
    integer :: i
    character(len=64), parameter, dimension(5) :: var_names = [ character(64) ::  &
        "temperature", "pressure", "density", "clv", "wind_u" ]

    config = atlas_Config()
    call config%set("type", "MultiFieldCreatorIFS");
    call config%set("ngptot", ngptot);
    call config%set("nproma", nproma);
    allocate(field_configs(size(var_names)))
    do i = 1, size(var_names)
        field_configs(i) = atlas_Config()
        call field_configs(i)%set("name", trim(var_names(i)))
    end do
    call field_configs(4)%set("nvar", 4) ! clv has four subvariables
    call config%set("fields", field_configs)

    call config%set("nlev", 0); ! surface fields
    call config%set("datatype", "real32");
    mfield_1 = atlas_MultiField(config)

    call config%set("nlev", 4); ! fields are 3d
    call config%set("datatype", "real64");
    mfield_2 = atlas_MultiField(config)

    fieldset_1 = mfield_1%fieldset()
    FCTEST_CHECK_EQUAL(mfield_1%size(), 5)
    FCTEST_CHECK_EQUAL(fieldset_1%size(), 5)

    fieldset_2 = atlas_FieldSet()
    call fieldset_2%add(mfield_1%fieldset())
    field = fieldset_2%field("density")
    call field%data(fdata_f2d)
    fdata_f2d(1,1) = 2.
    call field%rename("dens")

    ! check data access directly though multifield
    call fieldset_1%data("dens", fdata_f2d)
    fdata_f2d(1,1) = 3.

    ! check access to the renamed variable
    field = fieldset_1%field("dens")
    call field%data(fdata_f2d)
    FCTEST_CHECK_EQUAL(fdata_f2d(1,1), 3._c_float)

    ! check dimesionality
    fieldset_2 = mfield_2%fieldset()
    call fieldset_2%data("density", fdata_d3d)
    fdata_d3d(1,1,1) = 4.
    fieldset_2 = atlas_FieldSet()
    call fieldset_2%add(mfield_2%fieldset())
    field = fieldset_2%field("density")
    call field%data(fdata_d3d)
    FCTEST_CHECK_EQUAL(fdata_d3d(1,1,1), 4._c_double)
END_TEST


TEST( test_device_strides_on_cpu )
    implicit none

    type(atlas_MultiField)  :: mfield
    type(atlas_FieldSet)    :: fieldset
    type(atlas_Field)       :: field
    real(c_double), pointer :: fdata_d3d(:,:,:), fdata_d3d_ddata(:,:,:)

    integer, parameter :: nproma  = 6;
    integer, parameter :: nlev    = 5;
    integer, parameter :: ngptot  = 18;
    integer, parameter :: nblk    = (ngptot + nproma - 1) / nproma
    integer :: i, j, k

    character(len=64), parameter, dimension(2) :: var_names_dp = [ character(64) :: &
        "clv", "wind_u" ]

    mfield = atlas_MultiField(atlas_real(c_double), [nproma, nlev, -1, nblk], var_names_dp)

    FCTEST_CHECK_EQUAL(mfield%size(), 2)

    fieldset = mfield%fieldset()
    field = fieldset%field("clv")
    call field%data(fdata_d3d)
    do i = 1, field%shape(1)
      do j = 1, field%shape(2)
        do k = 1, field%shape(3)
          fdata_d3d(i,j,k) = real(1000*i + 100*j + 10*1 + k, c_double)
        end do
      end do
    end do
    call field%update_device()
    if (pluto%devices() == 0) then
      call field%device_data(fdata_d3d_ddata) ! We expect this to be host-resident, so we can read it
      do i = 1, field%shape(1)
        do j = 1, field%shape(2)
          do k = 1, field%shape(3)
            FCTEST_CHECK_EQUAL(fdata_d3d_ddata(i,j,k), fdata_d3d(i,j,k))
            fdata_d3d_ddata(i,j,k) = fdata_d3d_ddata(i,j,k) + 1_c_double
          end do
        end do
      end do
      call field%update_host()
      do i = 1, field%shape(1)
        do j = 1, field%shape(2)
          do k = 1, field%shape(3)
            FCTEST_CHECK_EQUAL(fdata_d3d(i,j,k), real(1000*i + 100*j + 10*1 + k + 1, c_double))
          end do
        end do
      end do
    endif
END_TEST


TEST( test_multifield_array_direct_constructor )
    implicit none

    type(atlas_MultiField)  :: mfield_1, mfield_2
    type(atlas_FieldSet)    :: fieldset_1, fieldset_2, fieldset_3
    type(atlas_Field)       :: field
    type(atlas_config)      :: config
    real(c_float), pointer  :: fdata_f2d(:,:)
    real(c_double), pointer :: fdata_d3d(:,:,:)

    integer, parameter :: nproma  = 6;
    integer, parameter :: nlev    = 5;
    integer, parameter :: ngptot  = 18;
    integer, parameter :: nblk    = (ngptot + nproma - 1) / nproma
    integer :: i, j, k
    character(len=64), parameter, dimension(3) :: var_names_sp = [ character(64) :: &
        "temperature ", "pressure", "density" ]
    character(len=64), parameter, dimension(2) :: var_names_dp = [ character(64) :: &
        "clv", "wind_u" ]

    mfield_1 = atlas_MultiField(atlas_real(c_float), [nproma, -1, nblk], var_names_sp)
    mfield_2 = atlas_MultiField(atlas_real(c_double), [nproma, nlev, -1, nblk], var_names_dp)

    FCTEST_CHECK_EQUAL(mfield_1%size(), 3)
    FCTEST_CHECK_EQUAL(mfield_2%size(), 2)

    fieldset_1 = mfield_1%fieldset()
    FCTEST_CHECK_EQUAL(fieldset_1%size(), 3)

    call fieldset_1%data("density", fdata_f2d)
    fdata_f2d(1,1) = 3._c_float

    fieldset_2 = atlas_FieldSet()
    call fieldset_2%add(mfield_2%fieldset())
    call fieldset_2%data("wind_u", fdata_d3d)
    fdata_d3d(1,1,1) = 4._c_double

    fieldset_3 = atlas_FieldSet()
    call fieldset_3%add(mfield_1%fieldset())
    call fieldset_3%add(mfield_2%fieldset())

    call fieldset_3%data("density", fdata_f2d)
    call fieldset_3%data("wind_u", fdata_d3d)
    FCTEST_CHECK_EQUAL(fdata_f2d(1,1), 3._c_float)
    FCTEST_CHECK_EQUAL(fdata_d3d(1,1,1), 4._c_double)

END_TEST


TEST( test_multifield_array_config_constuctor )
    implicit none

    type(atlas_MultiField)  :: mfield_1, mfield_2
    type(atlas_FieldSet)    :: fieldset_1, fieldset_2
    type(atlas_Field)       :: field
    type(atlas_config)      :: config
    integer, pointer :: fdata_int_2d(:,:)
    real(c_float), pointer  :: fdata_f2d(:,:), fdata_f3d(:,:,:)
    real(c_double), pointer :: fdata_d3d(:,:,:)

    integer, parameter :: nproma  = 16;
    integer, parameter :: nlev    = 1;
    integer, parameter :: nblk    = 200;
    type(atlas_Config), allocatable :: field_configs(:)
    integer :: i
    character(len=64), parameter, dimension(5) :: var_names = [ character(64) :: &
        "temperature", "pressure", "density", "clv", "wind_u" ]

    config = atlas_Config()
    call config%set("type", "MultiFieldCreatorArray");
    allocate(field_configs(size(var_names)))
    do i = 1, size(var_names)
        field_configs(i) = atlas_Config()
        call field_configs(i)%set("name", trim(var_names(i)))
    end do
    call field_configs(4)%set("nvar", 5) ! clv has four subvariables
    call config%set("fields", field_configs)

    ! surface fields
    call config%set("shape", [nproma, -1, nblk]);
    call config%set("datatype", "real32");
    mfield_1 = atlas_MultiField(config)

    ! fields are 3d
    call config%set("shape", [nproma, nlev, -1, nblk]);
    call config%set("datatype", "real64");
    mfield_2 = atlas_MultiField(config)

    fieldset_1 = mfield_1%fieldset()
    FCTEST_CHECK_EQUAL(mfield_1%size(), 9)
    FCTEST_CHECK_EQUAL(fieldset_1%size(), 9)

    fieldset_2 = atlas_FieldSet()
    call fieldset_2%add(mfield_1%fieldset())
    field = fieldset_2%field("density")
    call field%data(fdata_f2d)
    fdata_f2d(1,1) = 2.
    call field%rename("dens")

    ! check data access directly though multifield
    call fieldset_1%data("dens", fdata_f2d)
    fdata_f2d(1,1) = 3.

    ! check access to the renamed variable
    field = fieldset_1%field("dens")
    call field%data(fdata_f2d)
    FCTEST_CHECK_EQUAL(fdata_f2d(1,1), 3._c_float)

    ! check dimesionality
    fieldset_2 = mfield_2%fieldset()
    call fieldset_2%data("density", fdata_d3d)
    fdata_d3d(1,1,1) = 4.
    fieldset_2 = atlas_FieldSet()
    call fieldset_2%add(mfield_2%fieldset())
    field = fieldset_2%field("density")
    call field%data(fdata_d3d)
    FCTEST_CHECK_EQUAL(fdata_d3d(1,1,1), 4._c_double)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE
