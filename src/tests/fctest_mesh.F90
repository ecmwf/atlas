! (C) Copyright 1996-2015 ECMWF.
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
type(atlas_mesh_Nodes) :: nodes


type, extends(atlas_Config) :: atlas_FieldParametrisation
endtype

interface atlas_FieldParametrisation
  module procedure atlas_FieldConfig__ctor
end interface

contains

function atlas_FieldConfig__ctor(creator,ngptot,nproma,nlev,nvar,kind,datatype,shape,grid) result(params)
  use atlas_Config_c_binding
  type(atlas_FieldParametrisation) :: params
  character(len=*), optional, intent(in) :: creator
  integer, optional, intent(in) :: ngptot
  integer, optional, intent(in) :: nproma
  integer, optional, intent(in) :: nlev
  integer, optional, intent(in) :: nvar
  type(atlas_ReducedGrid), optional, intent(in) :: grid
  integer, optional, intent(in) :: kind
  character(len=*), optional, intent(in) :: datatype
  integer, optional, intent(in) :: shape(:)
  call params%reset_c_ptr( atlas__Config__new() )
  if( present(creator)   ) call params%set("creator"   ,creator  )
  if( present(ngptot)    ) call params%set("ngptot"    ,ngptot   )
  if( present(nproma)    ) call params%set("nproma"    ,nproma   )
  if( present(nlev)      ) call params%set("nlev"      ,nlev     )
  if( present(nvar)      ) call params%set("nvar"      ,nvar     )
  if( present(kind)      ) call params%set("kind"      ,kind     )
  if( present(datatype)  ) call params%set("datatype"  ,datatype )
  if( present(shape)     ) call params%set("shape"     ,shape    )
  call params%set("fortran",.True.) ! Let know that parameters have fortran style
end function

end module fctest_atlas_Mesh_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Mesh,fctest_atlas_Mesh_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
  mesh = atlas_Mesh()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call mesh%final()
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_function_space )

  write(*,*) "test_function_space starting"

  nodes = mesh%create_nodes(5)
  FCTEST_CHECK_EQUAL( nodes%size() , 5  )
  FCTEST_CHECK( nodes%has_field("partition") )
  FCTEST_CHECK( nodes%has_field("remote_idx") )
  call nodes%resize(10)
  call atlas_log%info( nodes%str() )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_name )
  type(atlas_Field) :: field

  field = atlas_Field("field",atlas_real(c_double),(/nodes%size()/))
  call nodes%add( field )
  FCTEST_CHECK_EQUAL( field%name() , "field" )
  call nodes%remove_field("field")
  call field%final() ! memory leak if not finalized!
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_owners)
  type(atlas_Field) :: f
  type(atlas_Field) :: f2
  type(atlas_State) :: state
  type(atlas_FieldSet) :: fields
  write(0,*) "test_field_owners"
  f = atlas_Field("field_test_owners",atlas_real(c_double),(/10/) )
  FCTEST_CHECK_EQUAL( f%owners() , 1 )
  state = atlas_State()
  call state%add(f)
  FCTEST_CHECK_EQUAL( f%owners() , 2 )

  f2 = state%field("field_test_owners")
  FCTEST_CHECK_EQUAL( f%owners() , 3 )
  call f2%final()

  FCTEST_CHECK_EQUAL( f%owners() , 2 )

  call state%remove("field_test_owners")
  FCTEST_CHECK_EQUAL( f%owners() , 1 )
  fields = atlas_FieldSet("fields")
  call fields%add(f)
  FCTEST_CHECK_EQUAL( f%owners() , 2 )

  call fields%final()
  FCTEST_CHECK_EQUAL( f%owners() , 1 )

  call f%final() ! memory leak if not finalized!
  call state%final()
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
  type(atlas_Field) field

  write(*,*) "test_field_metadata starting"

  field = atlas_Field("field_prop",atlas_real(c_float),(/1,nodes%size()/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
  call nodes%add(field)
  FCTEST_CHECK_EQUAL( field%owners() , 2 )
  metadata = field%metadata()
  call field%final()

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

  FCTEST_CHECK_EQUAL( int, 20 )
  FCTEST_CHECK_CLOSE( real32, real(0.1,kind=c_float), real(0.,kind=c_float) )
  FCTEST_CHECK_CLOSE( real64, real(0.2,kind=c_double), real(0.,kind=c_double) )
  FCTEST_CHECK_EQUAL( string, "hello world" )
  FCTEST_CHECK_EQUAL( arr_int32, (/1,2,3/) )
  FCTEST_CHECK_EQUAL( arr_real32, (/1.1_c_float,2.1_c_float,3.7_c_float/) )
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_size )
  integer, pointer :: fdata_int(:)
  real(c_float),  pointer :: fdata_real32(:)
  real(c_double), pointer :: fdata_real64(:)
  type(atlas_Field) :: field

  write(*,*) "test_field_size starting"

  field = atlas_Field("field_0",atlas_integer(),(/0,nodes%size()/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
  call nodes%add(field)
  FCTEST_CHECK_EQUAL( field%owners() , 2 )
  call field%data(fdata_int)
  FCTEST_CHECK_EQUAL( field%datatype() , "int32" )
  FCTEST_CHECK_EQUAL( size(fdata_int) , 0 )

  call field%final() ! Not necessary, following "=" will handle it
  write(0,*) "finalized field0"

  field = atlas_Field("field_1",atlas_real(c_float),(/1,nodes%size()/))
  call nodes%add(field)
  call field%data(fdata_real32)
  FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
  FCTEST_CHECK_EQUAL( size(fdata_real32) , 10 )

  call field%final() !Not necessary, following "=" will handle it

  field = atlas_Field("field_2",atlas_real(c_double),(/2,nodes%size()/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
  call nodes%add(field)
  FCTEST_CHECK_EQUAL( field%owners() , 2 )
  call field%data(fdata_real64)
  FCTEST_CHECK_EQUAL( field%name(), "field_2" )
  FCTEST_CHECK_EQUAL( field%datatype() , "real64" )
  FCTEST_CHECK_EQUAL( size(fdata_real64) , 20 )

  write(0,*) "Owners = ", field%owners()
  call field%attach()
  write(0,*) "Owners = ", field%owners()
  call field%attach()
  write(0,*) "Owners = ", field%owners()
  call field%detach()
  write(0,*) "Owners = ", field%owners()
  call field%detach()
  write(0,*) "Owners = ", field%owners()
  field = field
  write(0,*) "Owners = ", field%owners()
  call field%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_create_remove )
  real(c_double), pointer :: scalar(:)
  real(c_float), pointer :: vector(:,:)
  type(atlas_Field) :: field

  write(*,*) "test_create_remove starting"

  field = atlas_Field("bla",atlas_integer(),(/1,nodes%size()/))
  call nodes%add(field)
  FCTEST_CHECK_EQUAL( field%name(), "bla" )
  call field%final()

!  call nodes%remove_field("bla")

  field = atlas_Field("vector_field",atlas_real(c_float),(/3,nodes%size()/))
  call nodes%add(field)
  call field%data(vector)
  FCTEST_CHECK_EQUAL( size(vector),   30 )
  FCTEST_CHECK_EQUAL( size(vector,1), 3   )
  FCTEST_CHECK_EQUAL( size(vector,2), 10   )
!  call field%final()

  field = atlas_Field("scalar_field",atlas_real(c_double),(/1,nodes%size()/))
  call nodes%add(field)
  call field%data(scalar)
  FCTEST_CHECK_EQUAL( size(scalar),   10 )
  FCTEST_CHECK_EQUAL( size(scalar,1), 10  )
  FCTEST_CHECK_EQUAL( field%owners() , 2 )
  call field%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset )
  type(atlas_FieldSet) :: fieldset
  type(atlas_Field) :: afield
  type(atlas_Field) :: field

  write(*,*) "test_fieldset starting"

  fieldset = atlas_FieldSet("fieldset")

  afield = nodes%field("field_0")
  write(0,*) "field%owners() = ", afield%owners()
  call fieldset%add( nodes%field("field_0") )
  call fieldset%add( nodes%field("field_1") )
  call fieldset%add( nodes%field("field_2") )

  call fieldset%add( nodes%field("vector_field") )

  FCTEST_CHECK_EQUAL( fieldset%size(), 4 )

  field = fieldset%field(1)
  FCTEST_CHECK_EQUAL( field%name(), "field_0" )
  field = fieldset%field(2)
  FCTEST_CHECK_EQUAL( field%name(), "field_1" )
  field = fieldset%field(3)
  FCTEST_CHECK_EQUAL( field%name(), "field_2" )
  field = fieldset%field(4)
  FCTEST_CHECK_EQUAL( field%name(), "vector_field" )
  call fieldset%final()
  write(0,*) "test_fieldset end"
END_TEST

TEST( test_meshgen )
  type(atlas_ReducedGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_FunctionSpace) :: edges
  type(atlas_Field) :: field
  type(atlas_functionspace_Nodes) :: functionspace_nodes
  type(atlas_HaloExchange) :: halo_exchange
  integer(c_int), pointer :: ridx(:)
  real(c_double), pointer :: arr1d(:), arr2d(:,:)
  integer :: i, nnodes, nghost

  write(*,*) "test_meshgen starting"

  grid = atlas_ReducedGrid("N24")
  mesh = atlas_generate_mesh(grid)
  nodes = mesh%nodes()

!  call atlas_generate_reduced_gaussian_grid(rgg,"T63")
  call atlas_build_parallel_fields(mesh)
  call atlas_build_periodic_boundaries(mesh)
  call atlas_build_halo(mesh,1)
  call atlas_build_edges(mesh)
  call atlas_build_pole_edges(mesh)
  call atlas_build_median_dual_mesh(mesh)

  nnodes = nodes%size()

  field = nodes%field("remote_idx")
  call field%data(ridx)
  nghost = 0
  do i=1,nnodes
    if( ridx(i) /= i ) nghost = nghost + 1
  enddo

  write(0,*) "nghost =",nghost


  call atlas_log%info( nodes%str() )

  field = nodes%field("dual_volumes")
  call field%data(arr1d)
  call field%final()

  functionspace_nodes = atlas_functionspace_Nodes(mesh,1)
  halo_exchange = functionspace_nodes%get_halo_exchange()
  call halo_exchange%execute(arr1d)

  edges = mesh%function_space("edges")
  field = edges%field("dual_normals")
  call field%data(arr2d)
  call field%final()

  call atlas_write_gmsh(mesh,"testf2.msh")

  call atlas_write_load_balance_report(mesh,"N24_loadbalance.dat")
END_TEST

TEST( test_griddistribution )
  type(atlas_ReducedGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_GridDistribution) :: griddistribution

  integer, allocatable :: part(:)
  integer :: jnode

  grid = atlas_ReducedGrid("O16")
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
  call griddistribution%final()

  call atlas_write_gmsh(mesh,"testf3.msh")

  deallocate(part)
END_TEST


TEST( test_parametrisation )
  type(atlas_Config) :: params
  integer :: value
  character(len=:), allocatable :: valuestr
  logical :: found
  params = atlas_Config()

  if( .not. params%get("notexisting",value) ) then
    !call atlas_abort("notexisting not found",atlas_code_location(__FILE__,__LINE__))
  endif

  call params%set("value3",3)
  if( params%get("value3",value) ) then
    write(atlas_log%msg,*) "value = ",value; call atlas_log%info()
  endif

  found = params%get("value4",value)

  write(atlas_log%msg,*) "value = ",value; call atlas_log%info()

  call params%set("valuesttr","hello world")

  allocate( character(len=30) :: valuestr )
  valuestr = "goodbye"
  found = params%get("valuestr",valuestr)
  write(atlas_log%msg,*) "valuestr = ",valuestr; call atlas_log%info()

  call params%final()
END_TEST


TEST( test_fieldcreation )
  use fctest_atlas_Mesh_fixture , only : FieldParams => atlas_FieldParametrisation
  type(atlas_ReducedGrid) :: grid
  type(atlas_Field) :: field
  type(atlas_FieldParametrisation) :: params

  field = atlas_Field(FieldParams(creator="ArraySpec",shape=[10,137,1,30]))
  write(0,*) field%name(), field%size()
  call field%final()

  grid = atlas_ReducedGrid("O80")
  params = atlas_FieldParametrisation(creator="IFS",nproma=1024,ngptot=grid%npts(),nlev=137,nvar=1,kind=4)
  field = atlas_Field(params)
  call params%final()

  write(0,*) field%name(), field%size(), field%shape(), field%datatype(), field%bytes()
  call field%final()

! Idea:
!   field = atlas_Field([ &
!      & atlas_Param("creator","ArraySpec"),&
!      & atlas_Param("shape",[10,137,1,30]),&
!      & atlas_Param("kind","real64") ])

  call grid%final()
END_TEST

TEST( test_fv )
      type(atlas_ReducedGrid) :: grid
      type(atlas_Mesh) :: mesh
      type(atlas_GridDistribution) :: griddistribution
      type(atlas_functionspace_Nodes) :: nodes_fs

      type(atlas_FunctionSpace) :: edges
      type(atlas_mesh_Nodes) :: nodes

      integer, allocatable :: nloen(:)
      integer, allocatable :: part(:)
      integer :: halo_size = 1


      allocate(nloen(36))
      nloen(1:32) = 64

      ! Create a new Reduced Gaussian Grid based on a nloen array
      call atlas_log%info("Creating grid")
      grid = atlas_ReducedGaussianGrid( nloen(1:32) )

      ! Grid distribution: all points belong to partition 1
      allocate( part(grid%npts()) )
      part(:) = 1
      griddistribution = atlas_GridDistribution(part, part0=1)

      ! Generate mesh with given grid and distribution
      mesh = atlas_generate_mesh(grid,griddistribution)
      call griddistribution%final()

      ! Generate nodes function-space, with a given halo_size
      nodes_fs = atlas_functionspace_Nodes(mesh,halo_size)

      ! Create edge elements from surface elements
      call atlas_build_edges(mesh)
      call atlas_build_pole_edges(mesh)

      ! Generate edges function-space (This api will change soon)
      edges = mesh%function_space("edges")
      nodes = mesh%nodes()
      call atlas_build_edges_parallel_fields(edges,nodes)

      ! Build node to edge connectivity
      call atlas_build_node_to_edge_connectivity(mesh)

      ! Generate median-dual mesh, (dual_normals, dual_volumes)
      call atlas_build_median_dual_mesh(mesh)

      call mesh%final()
      call grid%final()
      call nodes_fs%final()

END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

