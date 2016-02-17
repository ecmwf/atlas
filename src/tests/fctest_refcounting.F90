! (C) Copyright 1996-2016 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.


#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_refcounting_fixture
use atlas_refcounted_module, only: atlas_RefCounted
use fctest
implicit none

type, extends(atlas_RefCounted) :: RefObj
contains
  procedure :: delete => RefObj__delete
  procedure :: copy => RefObj__copy
endtype

interface RefObj
  module procedure RefObj__constructor
end interface

contains


function create_obj(id) result(obj)
  type(RefObj) :: obj
  integer :: id
  obj = RefObj(id)
  call obj%return()
end function

function RefObj__constructor(id) result(obj)
  use atlas_mesh_c_binding
  type(RefObj) :: obj
  integer :: id
  write(0,*) "constructing obj ", id
  call obj%reset_c_ptr( atlas__Mesh__new() )
  !obj%id = id
end function

subroutine RefObj__delete(this)
  use atlas_mesh_c_binding
  class(RefObj), intent(inout) :: this
  write(0,*) "deleting obj"!,this%id
  if( .not. this%is_null() ) call atlas__Mesh__delete(this%c_ptr())
  call this%reset_c_ptr()
  write(0,*) "deleting obj"!,this%id, "done"
end subroutine


subroutine RefObj__copy(this,obj_in)
  class(RefObj), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
  write(0,*) "copy obj"
end subroutine

subroutine consume_new_obj(obj)
  type(RefObj) :: obj
  FCTEST_CHECK_EQUAL( obj%owners(), 0 )
  call consume_obj(obj)
end subroutine

subroutine consume_existing_obj(obj)
  type(RefObj) :: obj
  FCTEST_CHECK_EQUAL( obj%owners(), 1 )
  call consume_obj(obj)
end subroutine

subroutine consume_obj(obj)
  class(RefObj) :: obj
  call obj%attach()
  write(0,*) "Consume obj"!,obj%id
  call obj%final()
end subroutine

end module fctest_atlas_refcounting_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_refcounting,fctest_atlas_refcounting_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_ref )
  type(RefObj) :: obj, bjo
  obj = RefObj(1)
  FCTEST_CHECK_EQUAL( obj%owners(), 1 )
  !FCTEST_CHECK_EQUAL( obj%id, 1 )
  call consume_existing_obj(obj)
  call consume_new_obj(RefObj(2))

  !call obj%final() will be done in next statement upon assignment(=)
  obj = create_obj(3)
  FCTEST_CHECK_EQUAL( obj%owners(), 1 )
  bjo = obj
  FCTEST_CHECK_EQUAL( obj%owners(), 2 )
  obj = bjo
  FCTEST_CHECK_EQUAL( obj%owners(), 2 )
  call obj%final()
  call consume_existing_obj(bjo)
  call bjo%final()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

