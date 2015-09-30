! (C) Copyright 1996-2015 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.


#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_refcounting_fixture
use atlas_module
use iso_c_binding
use atlas_c_interop
use fctest
implicit none

type, extends(atlas_object) :: RefObj
  integer :: id
contains
  procedure :: finalize => RefObj__finalize
  procedure :: reset => RefObj__reset
  generic, public :: assignment(=) => reset
endtype

interface RefObj
  module procedure RefObj__constructor
end interface

contains

function create_obj(id) result(obj)
  type(RefObj) :: obj
  integer :: id
  obj = RefObj(id)
  call atlas_return(obj)
end function

function RefObj__constructor(id) result(obj)
  type(RefObj) :: obj
  integer :: id
  write(0,*) "constructing obj ", id
  call obj%reset_c_ptr( atlas__Mesh__new() )
  obj%id = id
end function

subroutine RefObj__finalize(obj)
  class(RefObj) :: obj
  if( obj%owners() > 0 ) call obj%detach()
  if( obj%owners() == 0 ) then
    write(0,*) "deleting obj",obj%id
    if( .not. obj%is_null() ) call atlas__Mesh__delete(obj%c_ptr())
    write(0,*) "reset_c_ptr",obj%id
    call obj%reset_c_ptr()
    write(0,*) "deleting obj",obj%id, "done"
  endif
end subroutine

subroutine RefObj__reset(obj_out,obj_in)
  use atlas_C_interop
  class(RefObj), intent(inout) :: obj_out
  class(RefObj), intent(in) :: obj_in
  if( obj_out /= obj_in ) then
    if( .not. obj_out%is_null() ) then
      write(0,*) "Finalizing obj",obj_out%id
      call obj_out%finalize()
    endif
    call obj_out%reset_c_ptr( obj_in%c_ptr() )
    if( .not. obj_out%is_null() ) call obj_out%attach()
    obj_out%id = obj_in%id
  endif
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
  write(0,*) "Consume obj",obj%id
  call obj%finalize()
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
  FCTEST_CHECK_EQUAL( obj%id, 1 )
  call consume_existing_obj(obj)
  call consume_new_obj(RefObj(2))

  !call obj%finalize() will be done in next statement upon assignment(=)
  obj = create_obj(3)
  FCTEST_CHECK_EQUAL( obj%owners(), 1 )
  bjo = obj
  FCTEST_CHECK_EQUAL( obj%owners(), 2 )
  obj = bjo
  FCTEST_CHECK_EQUAL( obj%owners(), 2 )
  call obj%finalize()
  call consume_existing_obj(bjo)
  call bjo%finalize()
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

