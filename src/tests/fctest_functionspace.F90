! (C) Copyright 1996-2015 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_FunctionSpace_Fixture
use atlas_module
use iso_c_binding
implicit none

contains

end module fctest_atlas_FunctionSpace_Fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_FunctionSpace,fctest_atlas_FunctionSpace_Fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_nodes )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_NodesFunctionSpace) :: fs
type(atlas_Field) :: field, template
type(atlas_Nodes) :: nodes
integer :: halo_size, nb_nodes
halo_size = 1

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fs = atlas_NodesFunctionSpace(mesh,halo_size)
nodes = fs%nodes()
nb_nodes = fs%nb_nodes()
write(atlas_log%msg,*) "nb_nodes = ",nb_nodes; call atlas_log%info()

field = fs%create_field("",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("field",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("",atlas_real(c_float),(/2/))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),(/2,2/))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field("",template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()


field = fs%create_global_field("",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_global_field("field",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_global_field("",atlas_real(c_float),(/2/))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_global_field("field",atlas_integer(c_int),(/2,2/))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_global_field("",template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_global_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

call fs%final()
call mesh%final()
call grid%final()

END_TEST

! -----------------------------------------------------------------------------


TEST( test_nodescolumns )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_NodesFunctionSpace) :: fs
type(atlas_Field) :: field, template
integer :: halo_size, levels
halo_size = 1
levels = 10

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fs = atlas_NodesFunctionSpace(mesh,halo_size)

!levels = fs%nb_levels()
write(atlas_log%msg,*) "nb_levels = ",levels; call atlas_log%info()

field = fs%create_field("",atlas_real(c_float),levels)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("field",atlas_real(c_float),levels)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("",atlas_real(c_float),levels,[2])
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),levels,[2,2])
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field("",template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()


field = fs%create_global_field("",atlas_real(c_float),levels)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_global_field("field",atlas_real(c_float),levels)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_global_field("",atlas_real(c_float),levels,[2])
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_global_field("field",atlas_integer(c_int),levels,[2,2])
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_global_field("",template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_global_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

call fs%final()
call mesh%final()
call grid%final()

END_TEST

! -----------------------------------------------------------------------------

TEST( test_collectives )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_NodesFunctionSpace) :: fs2d
type(atlas_Field) :: field, global, scal
real(c_float), pointer :: scalvalues(:)
real(c_float), pointer :: values(:,:)
real(c_float), pointer :: values3d(:,:,:,:)
real(c_float) :: minimum, maximum, sum, oisum, mean, stddev
real(c_float), allocatable :: minimumv(:), maximumv(:), sumv(:), oisumv(:), meanv(:), stddevv(:)
integer :: halo_size, levels
integer(ATLAS_KIND_GIDX) :: glb_idx
integer(ATLAS_KIND_GIDX), allocatable :: glb_idxv (:)
halo_size = 1
levels = 10

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fs2d = atlas_NodesFunctionSpace(mesh,halo_size)

field  = fs2d%create_field("",atlas_real(c_float),[2])
global = fs2d%create_global_field("",field)
scal   = fs2d%create_field("",atlas_real(c_float))

write(atlas_log%msg,*) "field:  rank",field%rank(), " shape [",field%shape(), "] size ", field%size();  call atlas_log%info()
write(atlas_log%msg,*) "global: rank",global%rank()," shape [",global%shape(),"] size ", global%size(); call atlas_log%info()

call fs2d%gather(field,global)
call fs2d%halo_exchange(field)
call fs2d%scatter(global,field)

call field%access_data(values)
call scal%access_data(scalvalues)
values = 2.
scalvalues = 2.

call atlas_log%info(fs2d%checksum(field))

values = atlas_mpi_rank()
scalvalues = atlas_mpi_rank()


call fs2d%minimum(scal,minimum)
call fs2d%maximum(scal,maximum)
write(atlas_log%msg,*) "min = ",minimum, " max = ",maximum; call atlas_log%info()

call fs2d%minimum(field,minimumv)
call fs2d%maximum(field,maximumv)
write(atlas_log%msg,*) "minv = ",minimumv, " maxv = ",maximumv; call atlas_log%info()

call fs2d%minimum_and_location(scal,minimum,glb_idx)
write(atlas_log%msg,*) "min = ",minimum, " gidx = ",glb_idx; call atlas_log%info()
call fs2d%maximum_and_location(scal,maximum,glb_idx)
write(atlas_log%msg,*) "max = ",maximum, " gidx = ",glb_idx; call atlas_log%info()

call fs2d%minimum_and_location(field,minimumv,glb_idxv)
write(atlas_log%msg,*) "minv = ",minimumv, " gidxv = ",glb_idxv; call atlas_log%info()
call fs2d%maximum_and_location(field,maximumv,glb_idxv)
write(atlas_log%msg,*) "minv = ",maximum, " gidxv = ",glb_idxv; call atlas_log%info()

call fs2d%sum(scal,sum)
call fs2d%order_independent_sum(scal,oisum)
write(atlas_log%msg,*) "sum = ",sum, " oisum = ",oisum; call atlas_log%info()

call fs2d%mean(scal,mean)
call fs2d%mean(field,meanv)
write(atlas_log%msg,*) "mean = ",mean, "meanv = ", meanv; call atlas_log%info()

call fs2d%mean_and_standard_deviation(scal,mean,stddev)
call fs2d%mean_and_standard_deviation(field,meanv,stddevv)
write(atlas_log%msg,*) "mean = ",mean, "meanv = ", meanv; call atlas_log%info()
write(atlas_log%msg,*) "stddev = ",stddev, "stddevv = ", stddevv; call atlas_log%info()

call scal%final()
call field%final()
call global%final()


field  = fs2d%create_field("",atlas_real(c_float),levels,[2,3])
global = fs2d%create_global_field("",field)

write(atlas_log%msg,*) "field:  rank",field%rank(), " shape [",field%shape(), "] size ", field%size();  call atlas_log%info()
write(atlas_log%msg,*) "global: rank",global%rank()," shape [",global%shape(),"] size ", global%size(); call atlas_log%info()

call fs2d%gather(field,global)
call fs2d%halo_exchange(field)
call fs2d%scatter(global,field)

call field%access_data(values3d)
values3d = 2.

call atlas_log%info(fs2d%checksum(field))
call field%final()
call global%final()
call fs2d%final()

call mesh%final()
call grid%final()

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

