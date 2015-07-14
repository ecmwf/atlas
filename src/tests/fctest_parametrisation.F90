! (C) Copyright 1996-2014 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_parametrisation_fixture
use atlas_module
implicit none

end module fctest_atlas_parametrisation_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_parametrisation,fctest_atlas_parametrisation_fixture)

! -----------------------------------------------------------------------------
TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_parametrisation )

  type(atlas_Parametrisation) :: params
  type(atlas_Parametrisation) :: nested
  type(atlas_Parametrisation) :: list(2)
  logical :: found
  integer :: intval
  integer :: j

  type(atlas_Parametrisation) :: anested
  type(atlas_Parametrisation), allocatable :: alist(:)

  ! --------------------- SET ------------------

  ! Corresponding YAML:
  ! {
  !   nested:
  !     list:
  !       -
  !         l1: 21
  !         l2: 22
  !       -
  !         l1: 21
  !         l2: 22
  !     n1: 11
  !     n2: 12
  !   p1: 1
  !   p2: 2
  ! }

  params = atlas_Parametrisation()
  call params%set("p1",1)
  call params%set("p2",2)

  nested = atlas_Parametrisation()
  call nested%set("n1",11)
  call nested%set("n2",12)

  do j=1,2
    list(j) = atlas_Parametrisation()
    call list(j)%set("l1",21)
    call list(j)%set("l2",22)
  enddo
  call nested%set("list",list)
  call atlas_delete(list)

  call params%set("nested",nested)
  call atlas_delete(nested)

  ! --------------------- JSON ------------------

  call atlas_log%info("params = "//params%json())

  ! --------------------- GET ------------------

  found = params%get("p1",intval)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL( intval , 1 )

  found = params%get("nested",anested)
  FCTEST_CHECK( found )
  found = anested%get("n1",intval)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL(intval, 11)

  found = anested%get("list",alist)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL( size(alist), 2 )

  found = alist(1)%get("l1",intval)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL(intval, 21)
  found = alist(1)%get("l2",intval)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL(intval, 22)

  found = alist(2)%get("l1",intval)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL(intval, 21)
  found = alist(2)%get("l2",intval)
  FCTEST_CHECK( found )
  FCTEST_CHECK_EQUAL(intval, 22)

  call atlas_delete(alist)
  call atlas_delete(anested)


  ! ---------------------------------------------

  call atlas_delete(params)

END_TEST

TEST(test_parametrisation_json_string)
 type(atlas_Parametrisation) :: params
 type(atlas_Parametrisation), allocatable :: records(:)
 type(atlas_JSON) :: json
 character (len=:), allocatable :: name
 integer :: age
 integer :: jrec
 json=atlas_JSON('{"records":['//&
  &       '{"name":"Joe",   "age":30},'//&
  &       '{"name":"Alison","age":43}' //&
  &    ']}')

 params = atlas_Parametrisation(json)
 call atlas_log%info(params%json())
 if( params%get("records",records) ) then
   do jrec=1,size(records)
     if( .not. records(jrec)%get("name",name) ) call atlas_abort("name not found")
     if( .not. records(jrec)%get("age",age) )   call atlas_abort("age not found")
     write(atlas_log%msg,'(2A,I0,A)') name," is ",age," years old"; call atlas_log%info()
  enddo
   call atlas_delete(records)
 endif
 call atlas_delete(params)
END_TEST

TEST(test_parametrisation_json_file)
 type(atlas_Parametrisation) :: params
 type(atlas_Parametrisation), allocatable :: records(:)
 type(atlas_Parametrisation) :: location
 type(atlas_JSON) :: json
 character (len=:), allocatable :: name, company, street, city
 integer :: age
 integer :: jrec

 ! Write a json file
 OPEN (UNIT=9 , FILE="fctest_parametrisation.json", STATUS='REPLACE')
 write(9,'(A)') '{"location":{"city":"Reading","company":"ECMWF","street":"Shinfield Road"},'//&
 &'"records":[{"age":42,"name":"Anne"},{"age":36,"name":"Bob"}]}'
 CLOSE(9)

 params = atlas_Parametrisation( atlas_PathName("fctest_parametrisation.json") )
 call atlas_log%info("params = "//params%json())

 if( params%get("records",records) ) then
   do jrec=1,size(records)
     if( .not. records(jrec)%get("name",name) ) call atlas_abort("name not found")
     if( .not. records(jrec)%get("age",age) )   call atlas_abort("age not found")
     write(atlas_log%msg,'(2A,I0,A)') name," is ",age," years old"; call atlas_log%info()
   enddo
   call atlas_delete(records)
 endif
 if( params%get("location",location) ) then
   call atlas_log%info("location = "//location%json())

   if( location%get("company",company) ) then
     write(0,*) "company = ",company
   endif
   if( location%get("street",street) ) then
     write(0,*) "street = ",street
   endif
   if( location%get("city",city) ) then
     write(0,*) "city = ",city
   endif
   call atlas_delete(location)
 endif
 call atlas_delete(params)

END_TEST
! -----------------------------------------------------------------------------

END_TESTSUITE

