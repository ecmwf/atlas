! (C) Copyright 2013 ECMWF.
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
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_nodes )
#if 1
type(atlas_StructuredGrid) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_NodeColumns) :: fs
type(atlas_Field) :: field, template
type(atlas_mesh_Nodes) :: nodes
integer :: halo_size, nb_nodes
halo_size = 1

grid = atlas_StructuredGrid("N24")
meshgenerator = atlas_MeshGenerator()
mesh = meshgenerator%generate(grid)
call meshgenerator%final()
fs = atlas_functionspace_NodeColumns(mesh,halo_size)
nodes = fs%nodes()
nb_nodes = fs%nb_nodes()
write(msg,*) "nb_nodes = ",nb_nodes; call atlas_log%info(msg)

field = fs%create_field(atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(name="field",kind=atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(atlas_real(c_float),variables=2)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),variables=2*2)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(template,name="field")
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()


field = fs%create_field(atlas_real(c_float),global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(name="field",kind=atlas_real(c_float),global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(atlas_real(c_float),variables=2,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),variables=2*2,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(template,name="field",global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

call fs%final()
call mesh%final()
call grid%final()
#else
#warning test test_nodes disabled
#endif
END_TEST

! -----------------------------------------------------------------------------


TEST( test_nodescolumns )
#if 1
type(atlas_StructuredGrid) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_NodeColumns) :: fs
type(atlas_Field) :: field, template
integer :: halo_size, levels
halo_size = 1
levels = 10

grid = atlas_StructuredGrid("N24")
meshgenerator = atlas_MeshGenerator()
mesh = meshgenerator%generate(grid)
call meshgenerator%final()
fs = atlas_functionspace_NodeColumns(mesh,halo_size)

!levels = fs%nb_levels()
write(msg,*) "nb_levels = ",levels; call atlas_log%info(msg)

field = fs%create_field(atlas_real(c_float),levels=levels)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(name="field",kind=atlas_real(c_float),levels=levels)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(atlas_real(c_float),levels=levels,variables=2)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),levels=levels,variables=2*2)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(template,name="field")
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()


field = fs%create_field(atlas_real(c_float),levels=levels,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(name="field",kind=atlas_real(c_float),levels=levels,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(atlas_real(c_float),levels=levels,variables=2,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),levels=levels,variables=2*2,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template,global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(template,name="field",global=.True.)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

fs = atlas_functionspace_NodeColumns(mesh,levels=5)
field = fs%create_field(atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

FCTEST_CHECK_EQUAL( fs%owners(), 1 )
call fs%final()

FCTEST_CHECK_EQUAL( mesh%owners(), 1 )
call mesh%final()

FCTEST_CHECK_EQUAL( grid%owners(), 1 )
call grid%final()
#else
#warning test test_nodescolumns disabled
#endif
END_TEST

! -----------------------------------------------------------------------------

TEST( test_collectives )
#if 1
use fckit_mpi_module
type(atlas_StructuredGrid) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_NodeColumns) :: fs2d
type(atlas_Field) :: field, global, scal
type(atlas_Metadata) :: metadata
type(fckit_mpi_comm) :: mpi
real(c_float), pointer :: scalvalues(:)
real(c_float), pointer :: values(:,:)
real(c_float), pointer :: values3d(:,:,:)
real(c_float) :: minimum, maximum, sum, oisum, mean, stddev
real(c_float), allocatable :: minimumv(:), maximumv(:), meanv(:), stddevv(:)
integer :: halo_size, levels
integer(ATLAS_KIND_GIDX) :: glb_idx
integer(ATLAS_KIND_GIDX), allocatable :: glb_idxv (:)
integer(c_int) :: test_broadcast
mpi = fckit_mpi_comm()
halo_size = 1
levels = 10

grid = atlas_StructuredGrid("N24")
meshgenerator = atlas_MeshGenerator()
mesh = meshgenerator%generate(grid)
call meshgenerator%final()
fs2d = atlas_functionspace_NodeColumns(mesh,halo_size)

field  = fs2d%create_field(kind=atlas_real(c_float),variables=2)
global = fs2d%create_field(field,global=.True.)
scal   = fs2d%create_field(kind=atlas_real(c_float))

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


field  = fs2d%create_field(kind=atlas_real(c_float),levels=levels,variables=2*3)
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
#endif
END_TEST




! -----------------------------------------------------------------------------


TEST( test_edges )
#if 1
type(atlas_StructuredGrid) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_functionspace_EdgeColumns) :: fs
type(atlas_Field) :: field, template
type(atlas_mesh_Edges) :: edges
integer :: halo_size, nb_edges
type( atlas_trace ) :: trace
type( atlas_trace ) :: trace_a
type( atlas_trace ) :: trace_b

type( atlas_Output ) :: gmsh

trace = atlas_Trace("fctest_functionspace.F90",__LINE__,"test_edges")
halo_size = 0

grid = atlas_StructuredGrid("N24")
meshgenerator = atlas_MeshGenerator()
mesh = meshgenerator%generate(grid)

gmsh = atlas_Output_gmsh("test_edges.msh")
call gmsh%write(mesh)
call gmsh%final()

FCTEST_CHECK_EQUAL( mesh%owners(), 1 )
edges = mesh%edges()
FCTEST_CHECK_EQUAL( edges%owners(), 3 ) ! Mesh holds 2 references (facets == edges)

trace_a = atlas_Trace("fctest_functionspace.F90",__LINE__,"EdgeColumns, no-levels")
fs = atlas_functionspace_EdgeColumns(mesh)
call trace_a%final()

FCTEST_CHECK_EQUAL( mesh%owners(), 2 )
FCTEST_CHECK_EQUAL( edges%owners(), 3 )
edges = fs%edges()
FCTEST_CHECK_EQUAL( edges%owners(), 3 )
nb_edges = fs%nb_edges()
write(msg,*) "nb_edges = ",nb_edges; call atlas_log%info(msg)
field = fs%create_field(atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(name="field",kind=atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(atlas_real(c_float),variables=2)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),variables=2*2)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(name="field",template=template)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

field = fs%create_field(atlas_real(c_float),levels=10,variables=2)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),variables=2*2,levels=10)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(name="field",template=template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()


field = fs%create_field(kind=atlas_real(c_float),global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(name="field", kind=atlas_real(c_float),global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()

field = fs%create_field(kind=atlas_real(c_float),variables=2,global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(name="field",kind=atlas_integer(c_int),variables=2*2,global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template,global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call field%final()

field = fs%create_field(template,name="field",global=.true.)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call field%final()
call template%final()

trace_b = atlas_Trace("fctest_functionspace.F90",__LINE__,"EdgeColumns, 5 levels")
fs = atlas_functionspace_EdgeColumns(mesh,levels=5)
call trace_b%final()
field = fs%create_field(atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )
call field%final()


call fs%final()
call edges%final()
call mesh%final()
call grid%final()

call trace%final()

#else
#warning test test_edges disabled
#endif
END_TEST

TEST( test_structuredcolumns )
#if 1
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_StructuredColumns) :: fs
type(atlas_functionspace) :: fs_base
integer(ATLAS_KIND_IDX) :: i, j
character(len=10) str

type(atlas_Field) :: field
type(atlas_Field) :: field_xy
type(atlas_Field) :: field_global_index
type(atlas_Field) :: field_index_j
real(8), pointer  :: xy(:,:), x(:)
integer(ATLAS_KIND_GIDX), pointer :: global_index(:)
integer(ATLAS_KIND_IDX), pointer  :: index_j(:)
integer, parameter :: XX=1

grid = atlas_StructuredGrid("O8")
fs = atlas_functionspace_StructuredColumns(grid,halo=2)
write(0,*) __LINE__

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

do j=fs%j_begin_halo(),fs%j_end_halo()
  write(str,'(I4,A)') j, ' : '
  call atlas_log%info(str,newl=.false.)
  do i=fs%i_begin_halo(j),fs%i_end_halo(j)
    write(str,'(I4)') global_index( fs%index(i,j) )
    call atlas_log%info(str,newl=.false.)
    x(fs%index(i,j)) = xy(XX,fs%index(i,j))
  enddo
  call atlas_log%info("",newl=.true.)
enddo

call fs%halo_exchange(field)

FCTEST_CHECK_EQUAL( field_xy%owners(), 2 )
fs = atlas_functionspace_StructuredColumns(grid,levels=5)

field = fs%create_field(atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )


write(0,*) "before: name = ", fs%name()
write(0,*) "before: owners = ", fs%owners()
fs_base = field%functionspace()
write(0,*) "after: name = " , fs%name()
write(0,*) "after: owners = " , fs%owners()

FCTEST_CHECK_EQUAL( field%owners(), 1 )
call field%final()
FCTEST_CHECK_EQUAL( field_xy%owners(), 1 )
call field_xy%final()
call field_global_index%final()
call fs%final()
call fs_base%final()
call grid%final()
#else
#warning test test_structuredcolumns disabled
#endif
END_TEST

TEST( test_pointcloud )
#if 1
implicit none
type(atlas_StructuredGrid) :: grid
type(atlas_functionspace_PointCloud) :: fs
type(atlas_functionspace) :: fs_base
type(atlas_trace) :: trace
character(len=10) str

type(atlas_Field) :: field, field2
type(atlas_Field) :: field_lonlat
real(8), pointer  :: lonlat(:,:), x(:)

trace = atlas_Trace("fctest_functionspace.F90",__LINE__,"test_pointcloud")

grid = atlas_Grid("O8")
fs = atlas_functionspace_PointCloud(grid)

field = fs%create_field(name="field",kind=atlas_real(8))
FCTEST_CHECK_EQUAL( field%owners(), 1 )
FCTEST_CHECK_EQUAL( field%levels(), 0 )
field_lonlat = fs%lonlat()
FCTEST_CHECK_EQUAL( field_lonlat%owners(), 2 )
call field%data(x)
call field_lonlat%data(lonlat)

FCTEST_CHECK_EQUAL( field_lonlat%owners(), 2 )
fs = atlas_functionspace_PointCloud(grid)

field = fs%create_field(atlas_real(c_float), levels=5)
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%kind() , atlas_real(c_float) )


write(0,*) "before: name = ", fs%name()
write(0,*) "before: owners = ", fs%owners()
fs_base = field%functionspace()
write(0,*) "after: name = " , fs%name()
write(0,*) "after: owners = " , fs%owners()

field2 = fs%create_field(name="field2",kind=atlas_real(8),levels=3,variables=2)

FCTEST_CHECK_EQUAL( field%shape(), ([5,int(grid%size())]) )
FCTEST_CHECK_EQUAL( field2%shape(), ([2,3,int(grid%size())]) )

#ifndef _CRAYFTN
FCTEST_CHECK_EQUAL( field%owners(), 1 )
#endif
call field%final()
#ifndef _CRAYFTN
FCTEST_CHECK_EQUAL( field_lonlat%owners(), 1 )
#endif
call field_lonlat%final()
call fs%final()
call fs_base%final()
call grid%final()
#else
#warning test test_pointcloud disabled
#endif
END_TEST

TEST( test_pointcloud_partition_remote )
#if 1
use fckit_module
implicit none

type(atlas_functionspace_PointCloud) :: fs
type(atlas_Field) :: fld_points
type(atlas_Field) :: fld_ghost
type(atlas_Field) :: fld_partition
type(atlas_Field) :: fld_remote_index
type(atlas_Field) :: fld_values
type(atlas_Field) :: fld_values_save
type(atlas_FieldSet) :: fset

type(atlas_functionspace) :: fs_base
type(atlas_trace) :: trace

integer i, k
integer rank
real space
real(c_double), allocatable :: point_values(:,:)
integer(c_int), dimension(4) :: ghost_values
integer(c_int), dimension(4) :: partition_index
integer(ATLAS_KIND_IDX), dimension(4) :: remote_index
real(c_double), pointer :: field_values(:,:)
real(c_double), pointer :: field_values_save(:,:)
real(c_double), pointer :: field_values_new(:,:)
real(c_double), parameter :: tol = 1.e-12_dp

trace = atlas_Trace("fctest_functionspace.F90",__LINE__,"test_pointcloud_partition_remote")
rank = fckit_mpi%rank()

! the idea here is that we have a latitudinal circle around the equator
! that is partitioned equally
!
! test here is for 4 PEs (but could be more!)
!   PE 0          PE 1        PE 2      PE 3
!
! Owned values
!    0  45       90 135      180 225   270 315
!
! Each PE has owned data followed by ghost data
! Ghost data is 1 before followed by 1 point after owned data ie.
!
! Correct halo_exchange for PE 0
! 0  45 | 315  90


space = 360.0 /(2.0 * fckit_mpi%size())
fset = atlas_FieldSet()

allocate(point_values(2, 4))

point_values = reshape((/space * (2 * rank), 0.0, &
                         space * (2 * rank + 1), 0.0, &
                         space * (2 * rank - 1), 0.0, &
                         space * (2 * rank + 2), 0.0/), shape(point_values))

if (point_values(1, 3) < 0.0) point_values(1, 3) = 360.0 + point_values(1, 3)
if (point_values(1, 4) > 360.0 - tol) point_values(1, 4) = 0.0d0
fld_points = atlas_Field("lonlat", point_values(:, :))
call fset%add(fld_points)

ghost_values = (/ 0, 0, 1, 1 /)
fld_ghost = atlas_Field("ghost", ghost_values(:))
call fset%add(fld_ghost)

partition_index = (/ rank, rank, rank - 1,  rank + 1 /)
if (partition_index(3) < 0) partition_index(3) = fckit_mpi%size() - 1
if (partition_index(4) == fckit_mpi%size()) partition_index(4) = 0
fld_partition = atlas_Field("partition", partition_index(:))
call fset%add(fld_partition)

remote_index = (/ 1, 2, 2, 1 /)
fld_remote_index = atlas_Field("remote_index", remote_index(:))
call fset%add(fld_remote_index)

fs = atlas_functionspace_PointCloud(fset)
fld_values = fs%create_field(name="values", kind=atlas_real(c_double), levels=1)
call fld_values%data(field_values)
fld_values_save = fs%create_field(name="values_save", kind=atlas_real(c_double), levels=1)
call fld_values_save%data(field_values_save)

do i = 1, 4
  do k = 1, 1
    field_values(k, i) = point_values(1, i)
    field_values_save(k, i) = point_values(1, i)
  end do
end do

call fs%halo_exchange(fld_values)
call fld_values%data(field_values_new)

do i = 1, 4
  do k = 1, 1
    FCTEST_CHECK_CLOSE(field_values_new(k, i), field_values_save(k, i), tol)
  end do
end do

call fs%adjoint_halo_exchange(fld_values)

do i = 1, 4
  do k = 1, 1
    if (i < 3) then
      FCTEST_CHECK_CLOSE(field_values_new(k, i), 2.0 * field_values_save(k, i), tol)
    else
      FCTEST_CHECK_CLOSE(field_values_new(k, i), 0.0_dp, tol)
    end if
end do
end do

#else
#warning test_pointcloud_partition_remote disabled
#endif
END_TEST




! -----------------------------------------------------------------------------

END_TESTSUITE

