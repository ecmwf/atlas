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

module fctest_atlas_Mesh_fixture
use atlas_module
use atlas_grids_module
use iso_c_binding
implicit none
type(atlas_Mesh) :: mesh
type(atlas_FunctionSpace) :: func_space
type(atlas_Field) :: field


type, extends(atlas_Metadata) :: atlas_FieldParametrisation
endtype

interface atlas_FieldParametrisation
  module procedure atlas_FieldParametrisation__ctor
end interface

contains

function atlas_FieldParametrisation__ctor(creator,ngptot,nproma,nlev,nvar,kind,data_type,shape,grid) result(params)
  type(atlas_FieldParametrisation) :: params
  character(len=*), optional, intent(in) :: creator
  integer, optional, intent(in) :: ngptot
  integer, optional, intent(in) :: nproma
  integer, optional, intent(in) :: nlev
  integer, optional, intent(in) :: nvar
  type(atlas_ReducedGrid), optional, intent(in) :: grid
  integer, optional, intent(in) :: kind
  character(len=*), optional, intent(in) :: data_type
  integer, optional, intent(in) :: shape(:)
  params%cpp_object_ptr = atlas__Metadata__new()
  if( present(creator)   ) call params%set("creator"   ,creator  )
  if( present(ngptot)    ) call params%set("ngptot"    ,ngptot   )
  if( present(nproma)    ) call params%set("nproma"    ,nproma   )
  if( present(nlev)      ) call params%set("nlev"      ,nlev     )
  if( present(nvar)      ) call params%set("nvar"      ,nvar     )
  if( present(kind)      ) call params%set("kind"      ,kind     )
  if( present(data_type) ) call params%set("data_type" ,data_type)
  if( present(shape)     ) call params%set("shape"     ,shape    )
  if( present(grid)      ) call params%set("grid"      ,grid     )
  call params%set("fortran",.True.) ! Let know that parameters have fortran style
end function

end module fctest_atlas_Mesh_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Mesh,fctest_atlas_Mesh_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
  call atlas_grids_load()
  mesh = atlas_Mesh()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_delete( mesh )
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_function_space )

  write(*,*) "test_function_space starting"

  call mesh%create_function_space( "nodes", "P1-C", 10 )
  func_space = mesh%function_space("nodes")
  CHECK_EQUAL( func_space%name() , "nodes"  )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_name )
  call func_space%create_field("field",1,atlas_real(c_double))
  field = func_space%field("field")
  CHECK_EQUAL( field%name() , "field" )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_metadata )
  integer :: int
  logical :: true, false
  real(c_float) :: real32
  real(c_double) :: real64
  character(len=:), allocatable :: string
  integer(c_int), allocatable :: arr_int32(:)
  real(c_float), allocatable :: arr_real32(:)
  type(atlas_Metadata) metadata
  type(atlas_Mesh) meshref

  write(*,*) "test_field_metadata starting"

  call func_space%create_field("field_prop",1,atlas_real(c_float))
  field = func_space%field("field_prop")

  metadata = field%metadata()

  call metadata%set("true",.True.)
  call metadata%set("false",.False.)
  call metadata%set("int",20)
  call metadata%set("real32", real(0.1,kind=c_float)   )
  call metadata%set("real64", real(0.2,kind=c_double) )
  call metadata%set("string", "hello world")
  call metadata%set("arr_int32", (/1,2,3/))
  call metadata%set("arr_int64", (/1_c_long,2_c_long,3_c_long/))
  call metadata%set("arr_real32", (/1.1_c_float,2.1_c_float,3.7_c_float/))
  call metadata%set("arr_real64", (/1.1_c_double,2.1_c_double,3.7_c_double/))
  call metadata%set("mesh", mesh )

  call metadata%get("true",true)
  call metadata%get("false",false)
  call metadata%get("int",int)
  call metadata%get("real32",real32)
  call metadata%get("real64",real64)
  call metadata%get("string",string)
  call metadata%get("arr_int64",arr_int32)
  call metadata%get("arr_real64",arr_real32)
  call metadata%get("mesh",meshref)

  call metadata%print(atlas_log%channel_info)

  write(atlas_log%msg,*) metadata%json()
  call atlas_log%info()
  !write(0,*) metadata%json()

  CHECK( true  .eqv. .True.  )
  CHECK( false .eqv. .False. )

  CHECK_EQUAL( int, 20 )
  CHECK_CLOSE( real32, real(0.1,kind=c_float), real(0.,kind=c_float) )
  CHECK_CLOSE( real64, real(0.2,kind=c_double), real(0.,kind=c_double) )
  CHECK_EQUAL( string, "hello world" )
  CHECK_EQUAL( arr_int32, (/1,2,3/) )
  CHECK_EQUAL( arr_real32, (/1.1,2.1,3.7/) )

END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_size )
  integer, pointer :: fdata_int(:)
  real(c_float),  pointer :: fdata_real32(:)
  real(c_double), pointer :: fdata_real64(:)

  write(*,*) "test_field_size starting"

  call func_space%create_field("field_0",0,atlas_integer())
  field = func_space%field("field_0")
  call field%access_data(fdata_int)
  CHECK_EQUAL( field%data_type() , "int32" )
  CHECK_EQUAL( size(fdata_int) , 0 )

  call func_space%create_field("field_1",1,atlas_real(c_float))
  field = func_space%field("field_1")
  call field%access_data(fdata_real32)
  CHECK_EQUAL( field%data_type() , "real32" )
  CHECK_EQUAL( size(fdata_real32) , 10 )

  call func_space%create_field("field_2",2,atlas_real(c_double))
  field = func_space%field("field_2")
  call field%access_data(fdata_real64)
  CHECK_EQUAL( field%name(), "field_2" )
  CHECK_EQUAL( field%data_type() , "real64" )
  CHECK_EQUAL( size(fdata_real64) , 20 )

END_TEST

! -----------------------------------------------------------------------------

TEST( test_create_remove )
  real(c_double), pointer :: scalar(:)
  real(c_float), pointer :: vector(:,:)

  write(*,*) "test_create_remove starting"

  call func_space%create_field("bla",1,atlas_integer())
  field = func_space%field("bla")
  CHECK_EQUAL( field%name(), "bla" )

!  call func_space%remove_field("bla")

  call func_space%create_field("vector_field",3,atlas_real(c_float))
  field = func_space%field("vector_field")
  call field%access_data(vector)
  CHECK_EQUAL( size(vector),   30 )
  CHECK_EQUAL( size(vector,1), 3   )
  CHECK_EQUAL( size(vector,2), 10   )

  call func_space%create_field("scalar_field",1,atlas_real(c_double))
  field = func_space%field("scalar_field")
  scalar => field%data1()
  CHECK_EQUAL( size(scalar),   10 )
  CHECK_EQUAL( size(scalar,1), 10  )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset )
  type(atlas_FieldSet) :: fieldset
  type(atlas_Field), allocatable :: fields(:)

  write(*,*) "test_fieldset starting"

  fieldset = atlas_FieldSet("fieldset")
  func_space = mesh%function_space("nodes")

  call fieldset%add_field( func_space%field("field_0") )
  call fieldset%add_field( func_space%field("field_1") )
  call fieldset%add_field( func_space%field("field_2") )

  call fieldset%add_field( func_space%field("vector_field") )

  CHECK_EQUAL( fieldset%size(), 4 )

  field = fieldset%field(1)
  CHECK_EQUAL( field%name(), "field_0" )
  field = fieldset%field(2)
  CHECK_EQUAL( field%name(), "field_1" )
  field = fieldset%field(3)
  CHECK_EQUAL( field%name(), "field_2" )
  field = fieldset%field(4)
  CHECK_EQUAL( field%name(), "vector_field" )


  call fieldset%get_array(fields)
  CHECK_EQUAL( size(fields), 4 )
  CHECK_EQUAL( fields(1)%name(), "field_0" )
  CHECK_EQUAL( fields(2)%name(), "field_1" )
  CHECK_EQUAL( fields(3)%name(), "field_2" )
  CHECK_EQUAL( fields(4)%name(), "vector_field" )

  call atlas_delete(fieldset)
END_TEST

TEST( test_meshgen )
  type(atlas_ReducedGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_FunctionSpace) :: nodes, edges
  type(atlas_Field) :: field
  integer, pointer :: bounds(:)
  integer(c_int), pointer :: ridx(:)
  real(c_double), pointer :: arr(:,:)
  integer :: i, nnodes, nghost

  write(*,*) "test_meshgen starting"

  grid = atlas_ReducedGrid("rgg.N24")
  mesh = atlas_generate_mesh(grid)

!  call atlas_generate_reduced_gaussian_grid(rgg,"T63")
  call atlas_build_parallel_fields(mesh)
  call atlas_build_periodic_boundaries(mesh)
  call atlas_build_halo(mesh,1)
  call atlas_build_edges(mesh)
  call atlas_build_pole_edges(mesh)
  call atlas_build_median_dual_mesh(mesh)

  nodes = mesh%function_space("nodes")
  call nodes%parallelise()
  bounds => nodes%shape()

  nnodes = bounds(2)

  field = nodes%field("remote_idx")
  call field%access_data(ridx)
  nghost = 0
  do i=1,nnodes
    if( ridx(i) /= i ) nghost = nghost + 1
  enddo

  write(0,*) "nghost =",nghost

  field = nodes%field("dual_volumes")
  call field%access_data(arr)
  call nodes%halo_exchange(arr)

  edges = mesh%function_space("edges")
  field = edges%field("dual_normals")
  call field%access_data(arr)


  !write(0,*) stride(ridx,1)

  !write(0,*) stride(arr,1), stride(arr,2), stride(arr,3)

  call atlas_write_gmsh(mesh,"testf2.msh")

  call atlas_write_load_balance_report(mesh,"N24_loadbalance.dat")
END_TEST

TEST( test_griddistribution )
  type(atlas_ReducedGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_GridDistribution) :: griddistribution

  integer, allocatable :: part(:)
  integer :: jnode

  grid = atlas_ReducedGrid("oct.N16")
  !grid = atlas_ReducedGrid("ll.128x64")
  !grid = atlas_LonLatGrid(128,64)

  allocate( part(grid%npts()) )
  do jnode=1,grid%npts()/3
    part(jnode) = 1
  enddo
  do jnode=grid%npts()/3+1,grid%npts()
    part(jnode) = 1
  enddo

  griddistribution = atlas_GridDistribution(part, part0=1)

  mesh = atlas_generate_mesh(grid,griddistribution)

  call atlas_write_gmsh(mesh,"testf3.msh")

  deallocate(part)
END_TEST

TEST( test_fieldcreation )
  use fctest_atlas_Mesh_fixture , only : FieldParams => atlas_FieldParametrisation
  type(atlas_ReducedGrid) :: grid
  type(atlas_Field) :: field
  type(atlas_FieldParametrisation) :: params

  field = atlas_Field(FieldParams(creator="ArraySpec",shape=[10,137,1,30]))
  write(0,*) field%name(), field%size()
  call atlas_delete(field)
  
  grid = atlas_ReducedGrid("oct.N80")
  params = atlas_FieldParametrisation(creator="IFS",nproma=1024,grid=grid,nlev=137,nvar=1,kind=4)
  field = atlas_Field(params)
  call atlas_delete(params)
  
  write(0,*) field%name(), field%size(), field%shape(), field%data_type(), field%bytes()
  call atlas_delete(field)
  
! Idea:
!   field = atlas_Field([ &
!      & atlas_Param("creator","ArraySpec"),&
!      & atlas_Param("shape",[10,137,1,30]),&
!      & atlas_Param("kind","real64") ])
  
  call atlas_delete(grid)
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

