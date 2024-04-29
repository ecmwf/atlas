! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
!
! @author Willem Deconinck
! @author Slavko Brdar

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_BlockStructuredColumns_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg
contains

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_BlockStructuredColumns,fcta_BlockStructuredColumns_fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

TEST( test_blockstructuredcolumns )
use atlas_functionspace_blockstructuredcolumns_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_functionspace) :: fs_base
integer(ATLAS_KIND_IDX) :: i, j, jbegin, jblk, jrof
character(len=256) str

type(atlas_Field) :: field
type(atlas_Field) :: field_xy
type(atlas_Field) :: field_global_index
type(atlas_Field) :: field_index_j
real(8), pointer  :: xy(:,:), x(:,:)
integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
integer(ATLAS_KIND_IDX), pointer  :: index_j(:)
integer(ATLAS_KIND_IDX) :: nblks, blksize
integer, parameter :: XX=1

grid = atlas_StructuredGrid("O16")
fs = atlas_functionspace_BlockStructuredColumns(grid,halo=0,nproma=12)

field = fs%create_field(name="field",kind=atlas_real(8))
FCTEST_CHECK_EQUAL( field%owners(), 1 )
field_xy = fs%xy()
FCTEST_CHECK_EQUAL( field_xy%owners(), 2 )
call field%data(x)
call field_xy%data(xy)
field_global_index = fs%global_index()
call field_global_index%data(global_index)

field_index_j = fs%index_j()
call field_index_j%data(index_j)

nblks=fs%nblks()
do jblk=1,nblks
  write(str,'(A,I4,A)') 'block', jblk, ' : '
  call atlas_log%info(str,newl=.false.)
  jbegin = fs%block_begin(jblk)
  blksize = fs%block_size(jblk)
  do jrof=1,blksize
    write(str,'(I6)') global_index( jbegin+jrof-1 )
    call atlas_log%info(str,newl=.false.)
    x(jrof,jblk) = xy(XX, jbegin+jrof-1)
  enddo
  call atlas_log%info("",newl=.true.)
enddo

!call fs%halo_exchange(field)

FCTEST_CHECK_EQUAL( field_xy%owners(), 2 )
fs = atlas_functionspace_BlockStructuredColumns(grid,levels=5)

field = fs%create_field(atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )


write(str,*) "before: name = ", fs%name(); call atlas_log%info(str)
write(str,*) "before: owners = ", fs%owners(); call atlas_log%info(str)
FCTEST_CHECK_EQUAL( fs%name(), "BlockStructuredColumns" )
fs_base = field%functionspace(); call atlas_log%info(str)
write(str,*) "after: name = " , fs_base%name(); call atlas_log%info(str)
write(str,*) "after: owners = " , fs_base%owners(); call atlas_log%info(str)
FCTEST_CHECK_EQUAL( fs_base%name(), "BlockStructuredColumns" )

FCTEST_CHECK_EQUAL( field%owners(), 1 )
call field%final()
FCTEST_CHECK_EQUAL( field_xy%owners(), 1 )
call field_xy%final()
call field_global_index%final()
call fs%final()
call fs_base%final()
call grid%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_collectives )
use fckit_mpi_module
use atlas_functionspace_blockstructuredcolumns_module
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_BlockStructuredColumns) :: fs
type(atlas_Field) :: field, global, scal
type(atlas_Config) :: config
type(atlas_Metadata) :: metadata
type(fckit_mpi_comm) :: mpi
integer(c_int) :: test_broadcast, halo_size, levels, nproma, nvar
call atlas_log%info("test_collectives")
mpi = fckit_mpi_comm()
halo_size = 0
levels = 10
nproma = 12
nvar = 2
FCTEST_CHECK_EQUAL( halo_size, 0 )
grid = atlas_StructuredGrid("O16")
fs = atlas_functionspace_BlockStructuredColumns(grid,nproma=nproma,levels=levels,halo=halo_size)

! type c_int
config = atlas_Config()
field  = fs%create_field(kind=atlas_integer(c_int),variables=2,options=config)
global = fs%create_field(field,global=.True.)
scal   = fs%create_field(kind=atlas_integer(c_int))
call fs%gather(field,global)
metadata = global%metadata()
if( mpi%rank() == 0 ) then
  call metadata%set("test_broadcast",123)
endif
call fs%scatter(global,field)
metadata = field%metadata()
call metadata%get("test_broadcast",test_broadcast)
FCTEST_CHECK_EQUAL( test_broadcast, 123 )

! type c_long
field  = fs%create_field(kind=c_long,variables=nvar)
global = fs%create_field(field,global=.True.)
scal   = fs%create_field(kind=c_long)
call fs%gather(field,global)
metadata = global%metadata()
if( mpi%rank() == 0 ) then
  call metadata%set("test_broadcast",123)
endif
call fs%scatter(global,field)
metadata = field%metadata()
call metadata%get("test_broadcast",test_broadcast)
FCTEST_CHECK_EQUAL( test_broadcast, 123 )

! type c_float
field  = fs%create_field(kind=atlas_real(c_float),variables=nvar)
global = fs%create_field(field,global=.True.)
scal   = fs%create_field(kind=atlas_real(c_float))
write(msg,*) "field:  rank",field%rank(), "        size ", field%size(), "         shape [",field%shape(), "]"
call atlas_log%info(msg)
write(msg,*) "nblk = ",fs%nblks(); call atlas_log%info(msg)
write(msg,*) "global: rank",global%rank(),"        size ", global%size(), "         shape [",global%shape(),"]"
call atlas_log%info(msg)
call fs%gather(field,global)
!call fs%halo_exchange(field)
metadata = global%metadata()
if( mpi%rank() == 0 ) then
  call metadata%set("test_broadcast",123)
  FCTEST_CHECK_EQUAL(global%shape(3), grid%size() )
  FCTEST_CHECK_EQUAL(global%shape(2), levels)
  FCTEST_CHECK_EQUAL(global%shape(1), nvar )
endif
call fs%scatter(global,field)
metadata = field%metadata()
call metadata%get("test_broadcast",test_broadcast)
FCTEST_CHECK_EQUAL( test_broadcast, 123 )

! type c_double
field  = fs%create_field(kind=atlas_real(c_double),variables=2)
global = fs%create_field(field,global=.True.)
scal   = fs%create_field(kind=atlas_real(c_double))
call fs%gather(field,global)
metadata = global%metadata()
if( mpi%rank() == 0 ) then
  call metadata%set("test_broadcast",123)
endif
call fs%scatter(global,field)
metadata = field%metadata()
call metadata%get("test_broadcast",test_broadcast)
FCTEST_CHECK_EQUAL( test_broadcast, 123 )
END_TEST

END_TESTSUITE

