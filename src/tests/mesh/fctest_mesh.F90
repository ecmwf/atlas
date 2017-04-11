! (C) Copyright 1996-2017 ECMWF.
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

module fcta_Mesh_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_Mesh_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Mesh,fcta_Mesh_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------


TEST( test_mesh_nodes )
implicit none

  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  integer(c_int) :: nb_nodes

  write(*,*) "test_function_space starting"
  mesh = atlas_Mesh()
  nodes = mesh%nodes()
  nb_nodes = nodes%size()
  FCTEST_CHECK_EQUAL( nb_nodes, 0 )
  FCTEST_CHECK_EQUAL( nodes%size() , 0_c_size_t  )
  FCTEST_CHECK( nodes%has_field("partition") )
  FCTEST_CHECK( nodes%has_field("remote_idx") )
  call nodes%resize(10_c_size_t)
  nb_nodes = nodes%size()
  FCTEST_CHECK_EQUAL( nb_nodes, 10 )
  FCTEST_CHECK_EQUAL( nodes%size() , 10_c_size_t  )
  call atlas_log%info( nodes%str() )

  call mesh%final()
  call nodes%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_name )
implicit none

  type(atlas_Field) :: field

  field = atlas_Field("field",atlas_real(c_double),(/10/))
  FCTEST_CHECK_EQUAL( field%name() , "field" )
  call field%final() ! memory leak if not finalized!
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_owners)
implicit none

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

TEST( test_field_size )
implicit none

  integer, pointer :: fdata_int(:,:)
  real(c_float),  pointer :: fdata_real32(:,:)
  real(c_double), pointer :: fdata_real64(:,:)
  type(atlas_Field) :: field

  write(*,*) "test_field_size starting"

  field = atlas_Field("field_0",atlas_integer(),(/0,10/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
  call field%data(fdata_int)
  FCTEST_CHECK_EQUAL( field%datatype() , "int32" )
  FCTEST_CHECK_EQUAL( size(fdata_int) , 0 )

  call field%final() ! Not necessary, following "=" will handle it
  write(0,*) "finalized field0"

  field = atlas_Field("field_1",atlas_real(c_float),(/1,10/))
  call field%data(fdata_real32)
  FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
  FCTEST_CHECK_EQUAL( size(fdata_real32) , 10 )

  call field%final() !Not necessary, following "=" will handle it

  field = atlas_Field("field_2",atlas_real(c_double),(/2,10/))
  FCTEST_CHECK_EQUAL( field%owners() , 1 )
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

TEST( test_field_wrapdata)
implicit none

  real(c_double), allocatable :: existing_data(:,:,:)
  real(c_double), pointer :: data(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: N=20
  integer(c_int) :: j
  write(0,*) "test_field_wrapdata"
  allocate( existing_data(2,10,N) )

  write(0,*) "line ", __LINE__

  ! Work with fields from here
  field = atlas_Field("wrapped",existing_data)

  write(0,*) "line ", __LINE__

  FCTEST_CHECK_EQUAL( field%rank()   , 3  )
  FCTEST_CHECK_EQUAL( field%size()   , 2*10*N )
  FCTEST_CHECK_EQUAL( field%shape(1) , 2  )
  FCTEST_CHECK_EQUAL( field%shape(2) , 10 )
  FCTEST_CHECK_EQUAL( field%shape(3) , N  )

  write(0,*) "line ", __LINE__

  FCTEST_CHECK_EQUAL( field%owners() , 1  )

  write(0,*) "line ", __LINE__

  call field%data(data)

  write(0,*) "line ", __LINE__

  do j=1,N
    data(1,1,j) = j
  enddo
  
  write(0,*) "line ", __LINE__
  
  call field%final()
  ! ... until here
  write(0,*) "line ", __LINE__

  ! Existing data is not deleted
  do j=1,N
    FCTEST_CHECK_EQUAL( existing_data(1,1,j), real(j,c_double) )
  enddo
  
  write(0,*) "line ", __LINE__
  
END_TEST

! -----------------------------------------------------------------------------

TEST( test_field_wrapdataslice)
implicit none

  real(c_double), allocatable :: existing_data(:,:,:,:)
  real(c_double), pointer :: data(:,:,:)
  type(atlas_Field) :: field
  integer(c_int) :: i,j,k,l
  write(0,*) "test_field_wrapdataslice"
  allocate( existing_data(4,3,2,5) )

  do i=1,4
    do j=1,3
      do k=1,2
        do l=1,5
          existing_data(i,j,k,l) = 1000*i+100*j+10*k+l
        enddo
      enddo
    enddo
  enddo

  field = atlas_Field(existing_data(:,:,1,:))
  FCTEST_CHECK_EQUAL( field%rank()   , 3  )
  FCTEST_CHECK_EQUAL( field%size()   , 4*3*5 )
  FCTEST_CHECK_EQUAL( field%shape(1) , 4  )
  FCTEST_CHECK_EQUAL( field%shape(2) , 3  )
  FCTEST_CHECK_EQUAL( field%shape(3) , 5  )

  call field%data(data)

  write(0,*) "Shape of field = ",shape(data)

  k=1
  do i=1,4
    do j=1,3
      do l=1,5
        FCTEST_CHECK_EQUAL(data(i,j,l) , real(1000*i+100*j+10*k+l,c_double) )
      enddo
    enddo
  enddo

  call field%final()

END_TEST

! -----------------------------------------------------------------------------

TEST( test_fieldset )
implicit none

  type(atlas_FieldSet) :: fieldset
  type(atlas_Field) :: field

  write(*,*) "test_fieldset starting"

  fieldset = atlas_FieldSet()

  field = atlas_Field("field_0",atlas_integer(),(/0,10/))
  call fieldset%add( field )

  field = atlas_Field("field_1",atlas_integer(),(/1,10/))
  call fieldset%add( field )

  field = atlas_Field("field_2",atlas_integer(),(/2,10/))
  call fieldset%add( field )

  FCTEST_CHECK_EQUAL( fieldset%size(), 3_c_size_t )

  field = fieldset%field(1)
  FCTEST_CHECK_EQUAL( field%name(), "field_0" )
  field = fieldset%field(2)
  FCTEST_CHECK_EQUAL( field%name(), "field_1" )
  field = fieldset%field(3)
  FCTEST_CHECK_EQUAL( field%name(), "field_2" )
  call fieldset%final()
  write(0,*) "test_fieldset end"
END_TEST

! -----------------------------------------------------------------------------

TEST( test_fv )
implicit none

      type(atlas_StructuredGrid) :: grid
      type(atlas_Mesh) :: mesh
      type(atlas_MeshGenerator) :: meshgenerator
      type(atlas_GridDistribution) :: griddistribution
      type(atlas_functionspace_NodeColumns) :: nodes_fs

      type(atlas_mesh_Edges) :: edges
      type(atlas_mesh_Nodes) :: nodes

      integer, allocatable :: nloen(:)
      integer, allocatable :: part(:)
      integer :: halo_size = 1


      type(atlas_Connectivity) :: node_to_node
      type(atlas_Connectivity) :: node_to_edge

      allocate(nloen(36))
      nloen(1:32) = 64

      ! Create a new Reduced Gaussian Grid based on a nloen array
      call atlas_log%info("Creating grid")
      grid = atlas_ReducedGaussianGrid( 32, nloen(1:32) )

      ! Grid distribution: all points belong to partition 1
      allocate( part(grid%npts()) )
      part(:) = 1
      griddistribution = atlas_GridDistribution(part, part0=1)

      ! Generate mesh with given grid and distribution
      meshgenerator = atlas_meshgenerator_Structured()
      mesh = meshgenerator%generate(grid,griddistribution)
      call griddistribution%final()

      ! Generate nodes function-space, with a given halo_size
      nodes_fs = atlas_functionspace_NodeColumns(mesh,halo_size)

      ! Create edge elements from surface elements
      call atlas_build_edges(mesh)
      call atlas_build_pole_edges(mesh)

      ! Generate edges function-space (This api will change soon)
      edges = mesh%edges()
      nodes = mesh%nodes()
      call atlas_build_edges_parallel_fields(mesh)

      ! Build node to edge connectivity
      call atlas_build_node_to_edge_connectivity(mesh)

      ! Generate median-dual mesh, (dual_normals, dual_volumes)
      call atlas_build_median_dual_mesh(mesh)


      node_to_node = atlas_Connectivity("node")
      call nodes%add(node_to_node)
      call node_to_node%final()

      node_to_node = nodes%connectivity("node")
      FCTEST_CHECK_EQUAL( node_to_node%rows(), 0_c_size_t )
      FCTEST_CHECK_EQUAL( node_to_node%name(),"node")

      node_to_node = nodes%connectivity("node")
      node_to_edge = nodes%connectivity("edge")
      
      write(0,*) "mesh.footprint (bytes) = ", mesh%footprint()

      call node_to_node%final()
      call mesh%final()
      call grid%final()
      call nodes_fs%final()

END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

