! (C) Copyright 1996-2016 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_FunctionSpace_fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
character(len=1024) :: msg
contains

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fcta_FunctionSpace,fcta_FunctionSpace_fxt)

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
type(atlas_grid_Structured) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_NodeColumns) :: fs
type(atlas_Field) :: field, template
type(atlas_mesh_Nodes) :: nodes
integer :: halo_size, nb_nodes
halo_size = 1

grid = atlas_grid_Structured("N24")
meshgenerator = atlas_meshgenerator_Structured()
mesh = meshgenerator%generate(grid)
call meshgenerator%final()
fs = atlas_functionspace_NodeColumns(mesh,halo_size)
nodes = fs%nodes()
nb_nodes = fs%nb_nodes()
write(msg,*) "nb_nodes = ",nb_nodes; call atlas_log%info(msg)

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

field = fs%create_field("",atlas_real(c_float),[2])
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),[2,2])
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


field = fs%create_field("",atlas_real(c_float),global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("field",atlas_real(c_float),global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("",atlas_real(c_float),[2],global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),[2,2],global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field("",template,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",template,global=.True.)
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
type(atlas_grid_Structured) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_NodeColumns) :: fs
type(atlas_Field) :: field, template
integer :: halo_size, levels
halo_size = 1
levels = 10

grid = atlas_grid_Structured("N24")
meshgenerator = atlas_meshgenerator_Structured()
mesh = meshgenerator%generate(grid)
call meshgenerator%final()
fs = atlas_functionspace_NodeColumns(mesh,halo_size)

!levels = fs%nb_levels()
write(msg,*) "nb_levels = ",levels; call atlas_log%info(msg)

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


field = fs%create_field("",atlas_real(c_float),levels,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("field",atlas_real(c_float),levels,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("",atlas_real(c_float),levels,[2],global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),levels,[2,2],global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field("",template,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",template,global=.True.)
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
use fckit_mpi_module
type(atlas_grid_Structured) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_NodeColumns) :: fs2d
type(atlas_Field) :: field, global, scal
type(atlas_Metadata) :: metadata
type(fckit_mpi_comm) :: mpi
real(c_float), pointer :: scalvalues(:)
real(c_float), pointer :: values(:,:)
real(c_float), pointer :: values3d(:,:,:,:)
real(c_float) :: minimum, maximum, sum, oisum, mean, stddev
real(c_float), allocatable :: minimumv(:), maximumv(:), sumv(:), oisumv(:), meanv(:), stddevv(:)
integer :: halo_size, levels
integer(ATLAS_KIND_GIDX) :: glb_idx
integer(ATLAS_KIND_GIDX), allocatable :: glb_idxv (:)
integer(c_int) :: test_broadcast
mpi = fckit_mpi_comm()
halo_size = 1
levels = 10

grid = atlas_grid_Structured("N24")
meshgenerator = atlas_meshgenerator_Structured()
mesh = meshgenerator%generate(grid)
call meshgenerator%final()
fs2d = atlas_functionspace_NodeColumns(mesh,halo_size)

field  = fs2d%create_field(atlas_real(c_float),[2])
global = fs2d%create_field(field,global=.True.)
scal   = fs2d%create_field(atlas_real(c_float))

write(msg,*) "field:  rank",field%rank(), " shape [",field%shape(), "] size ", field%size();  call atlas_log%info(msg)
write(msg,*) "global: rank",global%rank()," shape [",global%shape(),"] size ", global%size(); call atlas_log%info(msg)

call fs2d%gather(field,global)
call fs2d%halo_exchange(field)

metadata = global%metadata()
if( mpi%rank() == 0 ) then
  call metadata%set("test_broadcast",123)
endif

call fs2d%scatter(global,field)
metadata = field%metadata()
call metadata%get("test_broadcast",test_broadcast)
FCTEST_CHECK_EQUAL( test_broadcast, 123 )
call field%data(values)
call scal%data(scalvalues)
values = 2.
scalvalues = 2.

call atlas_log%info(fs2d%checksum(field))

values = mpi%rank()
scalvalues = mpi%rank()


call fs2d%minimum(scal,minimum)
call fs2d%maximum(scal,maximum)
write(msg,*) "min = ",minimum, " max = ",maximum; call atlas_log%info(msg)

call fs2d%minimum(field,minimumv)
call fs2d%maximum(field,maximumv)
write(msg,*) "minv = ",minimumv, " maxv = ",maximumv; call atlas_log%info(msg)

call fs2d%minimum_and_location(scal,minimum,glb_idx)
write(msg,*) "min = ",minimum, " gidx = ",glb_idx; call atlas_log%info(msg)
call fs2d%maximum_and_location(scal,maximum,glb_idx)
write(msg,*) "max = ",maximum, " gidx = ",glb_idx; call atlas_log%info(msg)

call fs2d%minimum_and_location(field,minimumv,glb_idxv)
write(msg,*) "minv = ",minimumv, " gidxv = ",glb_idxv; call atlas_log%info(msg)
call fs2d%maximum_and_location(field,maximumv,glb_idxv)
write(msg,*) "minv = ",maximum, " gidxv = ",glb_idxv; call atlas_log%info(msg)

call fs2d%sum(scal,sum)
call fs2d%order_independent_sum(scal,oisum)
write(msg,*) "sum = ",sum, " oisum = ",oisum; call atlas_log%info(msg)

call fs2d%mean(scal,mean)
call fs2d%mean(field,meanv)
write(msg,*) "mean = ",mean, "meanv = ", meanv; call atlas_log%info(msg)

call fs2d%mean_and_standard_deviation(scal,mean,stddev)
call fs2d%mean_and_standard_deviation(field,meanv,stddevv)
write(msg,*) "mean = ",mean, "meanv = ", meanv; call atlas_log%info(msg)
write(msg,*) "stddev = ",stddev, "stddevv = ", stddevv; call atlas_log%info(msg)

call scal%final()
call field%final()
call global%final()


field  = fs2d%create_field(atlas_real(c_float),levels,[2,3])
global = fs2d%create_field(field,global=.True.,owner=mpi%size()-1)

write(msg,*) "field:  rank",field%rank(), " shape [",field%shape(), "] size ", field%size();  call atlas_log%info(msg)
write(msg,*) "global: rank",global%rank()," shape [",global%shape(),"] size ", global%size(); call atlas_log%info(msg)

call fs2d%gather(field,global)
call fs2d%halo_exchange(field)
call fs2d%scatter(global,field)

call field%data(values3d)
values3d = 2.

call atlas_log%info(fs2d%checksum(field))
call field%final()
call global%final()
call fs2d%final()

call mesh%final()
call grid%final()

END_TEST




! -----------------------------------------------------------------------------


TEST( test_edges )
type(atlas_grid_Structured) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_EdgeColumns) :: fs
type(atlas_Field) :: field, template
type(atlas_mesh_Edges) :: edges
integer :: halo_size, nb_edges
halo_size = 0

grid = atlas_grid_Structured("N24")
meshgenerator = atlas_meshgenerator_Structured()
mesh = meshgenerator%generate(grid)
FCTEST_CHECK_EQUAL( mesh%owners(), 1 )
edges = mesh%edges()
FCTEST_CHECK_EQUAL( edges%owners(), 3 ) ! Mesh holds 2 references (facets == edges)
fs = atlas_functionspace_EdgeColumns(mesh)
FCTEST_CHECK_EQUAL( mesh%owners(), 2 )
FCTEST_CHECK_EQUAL( edges%owners(), 3 )
edges = fs%edges()
FCTEST_CHECK_EQUAL( edges%owners(), 3 )
nb_edges = fs%nb_edges()
write(msg,*) "nb_edges = ",nb_edges; call atlas_log%info(msg)

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

field = fs%create_field("",atlas_real(c_float),[2])
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),[2,2])
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


field = fs%create_field(atlas_real(c_float),global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field("field",atlas_real(c_float),global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(atlas_real(c_float),[2],global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",atlas_integer(c_int),[2,2],global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template,global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field("field",template,global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

call fs%final()
call edges%final()
call mesh%final()
call grid%final()

END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

