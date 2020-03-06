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

module fctest_atlas_Config_fixture
use atlas_module
use fckit_module, only : fckit_exception
implicit none

end module fctest_atlas_Config_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Config,fctest_atlas_Config_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_parametrisation )

  type(atlas_Config) :: params
  type(atlas_Config) :: nested
  type(atlas_Config), allocatable :: list(:)
  logical :: found
  integer :: intval
  integer :: j

  type(atlas_Config) :: anested
  type(atlas_Config), allocatable :: alist(:)

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

  params = atlas_Config()
  call params%set("p1",1)
  call params%set("p2",2)

  nested = atlas_Config()
  call nested%set("n1",11)
  call nested%set("n2",12)

  allocate( list(2) )
  do j=1,2
    list(j) = atlas_Config()
    call list(j)%set("l1",21)
    call list(j)%set("l2",22)
  enddo

  call nested%set("list",list)
  do j=1,2
    call list(j)%final()
  enddo

  call params%set("nested",nested)
  call nested%final()

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

  do j=1,size(alist)
    call alist(j)%final()
  enddo

  call anested%final()

  ! ---------------------------------------------

  call params%final()

END_TEST

! -----------------------------------------------------------------------------

TEST(test_parametrisation_json_string)
 type(atlas_Config) :: params
 type(atlas_Config), allocatable :: records(:)
 type(atlas_JSON) :: json
 character (len=:), allocatable :: name
 character(len=1024) :: msg
 integer :: age
 integer :: jrec
 json=atlas_JSON('{"records":['//&
  &       '{"name":"Joe",   "age":30},'//&
  &       '{"name":"Alison","age":43}' //&
  &    ']}')

 params = atlas_Config(json)
 call atlas_log%info(params%json())
 if( params%get("records",records) ) then
   do jrec=1,size(records)
     if( .not. records(jrec)%get("name",name) ) call fckit_exception%abort("name not found")
     if( .not. records(jrec)%get("age",age) )   call fckit_exception%abort("age not found")
     write(msg,'(2A,I0,A)') name," is ",age," years old"; call atlas_log%info(msg)
  enddo
  do jrec=1,size(records)
    call records(jrec)%final()
  enddo
 endif
 FCTEST_CHECK_EQUAL( params%owners() , 1 )
 call params%final()
END_TEST

TEST(test_parametrisation_json_file)
 type(atlas_Config) :: params
 type(atlas_Config), allocatable :: records(:)
 type(atlas_Config) :: location
 character (len=:), allocatable :: name, company, street, city
 integer :: age
 integer :: jrec
 character(len=1024) :: msg
 type(atlas_PathName) :: json_file

 ! Write a json file
 OPEN (UNIT=9 , FILE="fctest_parametrisation.json", STATUS='REPLACE')
 write(9,'(A)') '{"location":{"city":"Reading","company":"ECMWF","street":"Shinfield Road"},'//&
 &'"records":[{"age":42,"name":"Anne"},{"age":36,"name":"Bob"}]}'
 CLOSE(9)

 json_file = atlas_PathName("fctest_parametrisation.json")
 params = atlas_Config( json_file ) 
! params = atlas_Config( atlas_PathName("fctest_parametrisation.json") )
!     --> Does not work for XL compiler TODO: make reproducer
 call atlas_log%info("params = "//params%json())

 if( params%get("records",records) ) then
   do jrec=1,size(records)
     if( .not. records(jrec)%get("name",name) ) call fckit_exception%abort("name not found")
     if( .not. records(jrec)%get("age",age) )   call fckit_exception%abort("age not found")
     write(msg,'(2A,I0,A)') name," is ",age," years old"; call atlas_log%info(msg)
   enddo
   do jrec=1,size(records)
     call records(jrec)%final()
   enddo
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
   call location%final()
 endif
 call params%final()

END_TEST

! -----------------------------------------------------------------------------

TEST(test_json_file)
 type(atlas_JSON) :: json
 type(atlas_Config) :: config
 type(atlas_PathName) :: json_file
 ! Write a json file
 OPEN (UNIT=9 , FILE="fctest_parametrisation.json", STATUS='REPLACE')
 write(9,'(A)') '{"location":{"city":"Reading","company":"ECMWF","street":"Shinfield Road"},'//&
 &'"records":[{"age":42,"name":"Anne"},{"age":36,"name":"Bob"}]}'
 CLOSE(9)

 json_file = atlas_PathName( "fctest_parametrisation.json" )
 json = atlas_JSON( json_file )

 ! json = atlas_JSON( atlas_PathName("fctest_parametrisation.json") )
 !      --> Does not work with XL compiler TODO: make reproducer

 call atlas_log%info("json = "//json%str())

 config = atlas_Config( json )

 call atlas_log%info("config = "//config%json())

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

