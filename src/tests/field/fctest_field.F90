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

! -----------------------------------------------------------------------------

module fcta_Field_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Field,fcta_Field_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_field_name )
implicit none

  type(atlas_Field) :: field

  field = atlas_Field("field",atlas_real(c_double),(/10/))
  FCTEST_CHECK_EQUAL( field%name() , "field" )
  call field%final() ! memory leak if not finalized!
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_owners)
implicit none

  type(atlas_Field) :: f
  type(atlas_Field) :: f2
  type(atlas_State) :: state
  type(atlas_FieldSet) :: fields
  write(0,*) "test_field_owners"
  f = atlas_Field("field_test_owners",atlas_real(c_double),(/10/) )
  FCTEST_CHECK_EQUAL( f%owners() , 1 )
  state = atlas_State()
  call state%add(f)
  FCTEST_CHECK_EQUAL( f%owners() , 2 )

  f2 = state%field("field_test_owners")
  FCTEST_CHECK_EQUAL( f%owners() , 3 )
  call f2%final()

  FCTEST_CHECK_EQUAL( f%owners() , 2 )

  call state%remove("field_test_owners")
  FCTEST_CHECK_EQUAL( f%owners() , 1 )
  fields = atlas_FieldSet("fields")

  FCTEST_CHECK_EQUAL( fields%name() , "fields" )

  call fields%add(f)
  FCTEST_CHECK_EQUAL( f%owners() , 2 )

  call fields%final()
  FCTEST_CHECK_EQUAL( f%owners() , 1 )

  call f%final() ! memory leak if not finalized!
  call state%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_size )
implicit none

  integer, pointer :: fdata_int(:,:)
  real(c_float),  pointer :: fdata_real32(:,:)
  real(c_double), pointer :: fdata_real64(:,:)
  type(atlas_Field) :: field

  write(*,*) "test_field_size starting"

  field = atlas_Field("field_0",atlas_integer(),(/0,10/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
  call field%data(fdata_int)
  FCTEST_CHECK_EQUAL( field%datatype() , "int32" )
  FCTEST_CHECK_EQUAL( size(fdata_int) , 0 )

  call field%final() ! Not necessary, following "=" will handle it
  write(0,*) "finalized field0"

  field = atlas_Field("field_1",atlas_real(c_float),(/1,10/))
  call field%data(fdata_real32)
  FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
  FCTEST_CHECK_EQUAL( size(fdata_real32) , 10 )

  call field%final() !Not necessary, following "=" will handle it

  field = atlas_Field("field_2",atlas_real(c_double),(/2,10/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
  call field%data(fdata_real64)
  FCTEST_CHECK_EQUAL( field%name(), "field_2" )
  FCTEST_CHECK_EQUAL( field%datatype() , "real64" )
  FCTEST_CHECK_EQUAL( size(fdata_real64) , 20 )

  write(0,*) "Owners = ", field%owners()
  call field%attach()
  write(0,*) "Owners = ", field%owners()
  call field%attach()
  write(0,*) "Owners = ", field%owners()
  call field%detach()
  write(0,*) "Owners = ", field%owners()
  call field%detach()
  write(0,*) "Owners = ", field%owners()
  field = field
  write(0,*) "Owners = ", field%owners()
  call field%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset )
implicit none

  type(atlas_FieldSet) :: fieldset, fieldset_2
  type(atlas_Field) :: field, field_2

  write(*,*) "test_fieldset starting"

  fieldset = atlas_FieldSet()

  field = atlas_Field("field_0",atlas_integer(),(/0,10/))
  call fieldset%add( field )

  field = atlas_Field("field_1",atlas_integer(),(/1,10/))
  call fieldset%add( field )

  field = atlas_Field("field_2",atlas_integer(),(/2,10/))
  call fieldset%add( field )

  FCTEST_CHECK_EQUAL( fieldset%size(), 3 )

  field = fieldset%field(1)
  FCTEST_CHECK_EQUAL( field%name(), "field_0" )
  field = fieldset%field(2)
  FCTEST_CHECK_EQUAL( field%name(), "field_1" )
  field = fieldset%field(3)
  FCTEST_CHECK_EQUAL( field%name(), "field_2" )

  fieldset_2 = atlas_FieldSet()
  call fieldset_2%add(fieldset)
  field_2 = fieldset_2%field("field_0")
  call field_2%rename("field_00")
  FCTEST_CHECK(fieldset%has("field_00"))
  field = fieldset%field("field_00")
  FCTEST_CHECK_EQUAL(field%name(), "field_00")

  call fieldset%final()
  write(0,*) "test_fieldset end"

END_TEST

TEST( test_field_aligned )
    implicit none
    type(atlas_Field) :: field
    real(c_double), pointer :: view(:,:,:)
    field = atlas_Field("field_3",atlas_real(c_double),(/3,5,10/),alignment=4)
    call field%data(view)
    FCTEST_CHECK_EQUAL( size(view,1) , 3 )
    FCTEST_CHECK_EQUAL( size(view,2) , 5 )
    FCTEST_CHECK_EQUAL( size(view,3) , 10 )
    FCTEST_CHECK_EQUAL( field%shape(1), 3 )
    FCTEST_CHECK_EQUAL( field%shape(2), 5 )
    FCTEST_CHECK_EQUAL( field%shape(3), 10 )
    FCTEST_CHECK_EQUAL( field%stride(1), 1 )
    FCTEST_CHECK_EQUAL( field%stride(2), 4 )
    FCTEST_CHECK_EQUAL( field%stride(3), 4*5 )
    call field%final()
END_TEST

TEST( test_fieldset_slice )
    implicit none
    type(atlas_FieldSet) :: fieldset
    type(atlas_Field) :: field
    integer, pointer :: view1d(:), view3d(:,:,:)
    integer, pointer :: slice0d, slice2d(:,:)

    ! slicing of a three-dimensional field
    field = atlas_Field("field_4",atlas_integer(),(/3,5,10/))
    call field%data(view3d)
    view3d(1,2,3) = 123
    call field%data(slice2d,3)
    FCTEST_CHECK_EQUAL( size(slice2d,1) , 3 )
    FCTEST_CHECK_EQUAL( size(slice2d,2) , 5 )
    FCTEST_CHECK_EQUAL( slice2d(1,2), 123 )
    slice2d(1,2) = slice2d(1,2) * 2
    call field%data(view3d)
    FCTEST_CHECK_EQUAL( view3d(1,2,3), 246 )

    ! slicing of a one-dimensional field
    field = atlas_Field("field_5",atlas_integer(),(/3/))
    call field%data(view1d)
    view1d(2) = 123
    call field%data(slice0d,2)
    FCTEST_CHECK_EQUAL( slice0d, 123 )
    slice0d = slice0d * 2
    call field%data(view1d)
    FCTEST_CHECK_EQUAL( view1d(2), 246 )
    call field%final()

    ! slicing of a three-dimensional field through a fieldset by name and by idx
    fieldset = atlas_FieldSet()
    field = atlas_Field("field_6",atlas_integer(),(/3,5,10/))
    call fieldset%add(field)
    field = atlas_Field("field_7",atlas_integer(),(/3/))
    call fieldset%add(field)
    call fieldset%data("field_6", view3d)
    view3d(1,2,3) = 122
    call fieldset%data(1, view3d)
    view3d(1,2,3) = 123
    call fieldset%data("field_6", slice2d, 3)
    slice2d(1,2) = slice2d(1,2) + 1
    call fieldset%data(1, slice2d, 3)
    slice2d(1,2) = slice2d(1,2) - 1
    FCTEST_CHECK_EQUAL( size(slice2d,1) , 3 )
    FCTEST_CHECK_EQUAL( size(slice2d,2) , 5 )
    FCTEST_CHECK_EQUAL( slice2d(1,2), 123 )
    slice2d(1,2) = slice2d(1,2) * 2
    call fieldset%data('field_6', view3d)
    FCTEST_CHECK_EQUAL( view3d(1,2,3), 246 )

    ! slicing of a one-dimensional field through a fieldset by name and by idx
    call fieldset%data('field_7', view1d)
    view1d(2) = 122
    call fieldset%data(2, view1d)
    view1d(2) = 123
    call field%data(slice0d, 2)
    FCTEST_CHECK_EQUAL( slice0d, 123 )
    slice0d = slice0d * 2
    call field%data(view1d)
    FCTEST_CHECK_EQUAL( view1d(2), 246 )
    call field%final()
    call fieldset%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

