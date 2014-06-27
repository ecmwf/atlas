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

module fctest_atlas_mesh_fixture
use atlas_module
use atlas_mpl_module
use iso_c_binding
implicit none
type(Mesh_type) :: mesh
type(FunctionSpace_type) :: func_space
type(Field_type) :: field
end module fctest_atlas_mesh_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_mesh,fctest_atlas_mesh_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call MPL_Init()
  mesh = new_Mesh()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call delete( mesh )
  call MPL_Finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_function_space )
  call mesh%add_function_space( new_FunctionSpace("nodes", "P1-C", 10) )
  func_space = mesh%function_space("nodes")
  CHECK_EQUAL( func_space%name() , "nodes"  )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_name )
  call func_space%create_field("field",1,real_kind(c_double))
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
  type(MetaData_type) metadata
  call func_space%create_field("field_prop",1,real_kind(c_float))
  field = func_space%field("field_prop")

  metadata = field%metadata()

  call metadata%add("true",.True.)
  call metadata%add("false",.False.)
  call metadata%add("int",20)
  call metadata%add("real32", real(0.1,kind=c_float)   )
  call metadata%add("real64", real(0.2,kind=c_double) )
  call metadata%add("string", "hello world")

  call metadata%get("true",true)
  call metadata%get("false",false)
  call metadata%get("int",int)
  call metadata%get("real32",real32)
  call metadata%get("real64",real64)
  call metadata%get("string",string)
  
  CHECK( true  .eqv. .True.  )
  CHECK( false .eqv. .False. )

  CHECK_EQUAL( int, 20 )
  CHECK_CLOSE( real32, real(0.1,kind=c_float), real(0.,kind=c_float) )
  CHECK_CLOSE( real64, real(0.2,kind=c_double), real(0.,kind=c_double) )
  CHECK_EQUAL( string, "hello world" )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_size )
  integer, pointer :: fdata_int(:)
  real(kind=c_float),  pointer :: fdata_real32(:)
  real(kind=c_double), pointer :: fdata_real64(:)
  call func_space%create_field("field_0",0,integer_kind())
  field = func_space%field("field_0")
  call field%access_data(fdata_int)
  CHECK_EQUAL( field%data_type() , "int32" )
  CHECK_EQUAL( size(fdata_int) , 0 )

  call func_space%create_field("field_1",1,real_kind(c_float))
  field = func_space%field("field_1")
  call field%access_data(fdata_real32)
  CHECK_EQUAL( field%data_type() , "real32" )
  CHECK_EQUAL( size(fdata_real32) , 10 )

  call func_space%create_field("field_2",2,real_kind(c_double))
  field = func_space%field("field_2")
  call field%access_data(fdata_real64)
  CHECK_EQUAL( field%name(), "field_2" )
  CHECK_EQUAL( field%data_type() , "real64" )
  CHECK_EQUAL( size(fdata_real64) , 20 )

END_TEST

! -----------------------------------------------------------------------------

TEST( test_create_remove )
  call func_space%create_field("bla",1,integer_kind())
  field = func_space%field("bla")
  CHECK_EQUAL( field%name(), "bla" )

  call func_space%remove_field("bla")
END_TEST

! -----------------------------------------------------------------------------

TEST( test_prismatic_function_space )
  real(kind=c_double), pointer :: scalar(:,:)
  real(kind=c_float), pointer :: vector(:,:,:)
  call mesh%add_function_space( new_PrismaticFunctionSpace("prismatic", "P1-C", 5, 10 ) ) 

  func_space = mesh%function_space("prismatic")
  call func_space%create_field("vector_field",3,real_kind(c_float))
  field = func_space%field("vector_field")
  call field%access_data(vector)
  CHECK_EQUAL( size(vector),   150 )
  CHECK_EQUAL( size(vector,1), 3   )
  CHECK_EQUAL( size(vector,2), 5   )
  CHECK_EQUAL( size(vector,3), 10  )

  call func_space%create_field("scalar_field",1,real_kind(c_double))
  field = func_space%field("scalar_field")
  scalar => field%data2()
  CHECK_EQUAL( size(scalar),   50 )
  CHECK_EQUAL( size(scalar,1), 5  )
  CHECK_EQUAL( size(scalar,2), 10 )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset )
  type(FieldSet_type) :: fieldset
  type(Field_type), allocatable :: fields(:)
  fieldset = new_FieldSet("fieldset")
  func_space = mesh%function_space("nodes")

  call fieldset%add_field( func_space%field("field_0") )
  call fieldset%add_field( func_space%field("field_1") )
  call fieldset%add_field( func_space%field("field_2") )

  func_space = mesh%function_space("prismatic")
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

  call delete(fieldset)
END_TEST

TEST( test_meshgen )

  type(Mesh_type) :: rgg
  type(FunctionSpace_type) :: nodes, edges
  type(Field_type) :: field
  integer, pointer :: bounds(:)
  integer(c_int), pointer :: ridx(:)
  real(c_double), pointer :: arr(:,:)
  integer :: i, nnodes, nghost
  call atlas_generate_latlon_grid(rgg,60,30)

!  call atlas_generate_reduced_gaussian_grid(rgg,"T63")
  call atlas_build_parallel_fields(rgg)
  call atlas_build_periodic_boundaries(rgg)
  call atlas_build_halo(rgg,1)
  call atlas_build_edges(rgg)
  call atlas_build_pole_edges(rgg)
  call atlas_build_dual_mesh(rgg)

  nodes = rgg%function_space("nodes")
  call nodes%parallelise()
  bounds => nodes%bounds()

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

  edges = rgg%function_space("edges")
  field = edges%field("dual_normals")
  call field%access_data(arr)


  write(0,*) stride(ridx,1)

  write(0,*) stride(arr,1), stride(arr,2), stride(arr,3)

  call atlas_write_gmsh(rgg,"testf2.msh")

  call atlas_write_load_balance_report(rgg,"T63_loadbalance.dat")
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

