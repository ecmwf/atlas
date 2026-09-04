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

!!! WARNING ! THIS IS DEPRECATED AND SHOULD NOT BE USED AS EXAMPLE !!!!

TEST( test_mesh_nodes )
implicit none

  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  integer(c_int) :: nb_nodes

  call atlas_log%info( "--- test_mesh_nodes" )
  mesh = atlas_Mesh()
  nodes = mesh%nodes()
  nb_nodes = nodes%size()
  FCTEST_CHECK_EQUAL( nb_nodes, 0 )
  FCTEST_CHECK_EQUAL( nodes%size() , 0  )
  FCTEST_CHECK( nodes%has_field("partition") )
  FCTEST_CHECK( nodes%has_field("remote_idx") )
  call nodes%resize(10)
  nb_nodes = nodes%size()
  FCTEST_CHECK_EQUAL( nb_nodes, 10 )
  FCTEST_CHECK_EQUAL( nodes%size() , 10  )
  call atlas_log%info( nodes%str() )

  call mesh%final()
  call nodes%final()
END_TEST

TEST( test_mesh_grid_distribution )
implicit none

  type(atlas_Grid) :: grid
  type(atlas_GridDistribution) :: distribution
  type(atlas_Partitioner) :: partitioner
  type(atlas_Mesh) :: mixed_mesh
  type(atlas_mesh_Nodes) :: mixed_nodes
  type(atlas_mesh_Cells) :: mixed_cells
  type(atlas_Elements) :: mixed_quads, mixed_triangles
  type(atlas_ElementType) :: element_type

  call atlas_log%info( "--- test_mesh_grid_distribution" )

  grid = atlas_Grid("N16")
  partitioner = atlas_Partitioner("serial")
  distribution = atlas_GridDistribution(grid, partitioner)

  mixed_mesh = atlas_Mesh(grid, distribution)

  mixed_nodes = mixed_mesh%nodes()
  FCTEST_CHECK( mixed_nodes%size() > 0 )

  mixed_cells = mixed_mesh%cells()
  FCTEST_CHECK_EQUAL( mixed_cells%nb_types(), 2 )
  mixed_quads = mixed_cells%elements(1)
  mixed_triangles = mixed_cells%elements(2)
  FCTEST_CHECK_EQUAL( mixed_nodes%size(), 1720 )
  FCTEST_CHECK_EQUAL( mixed_quads%size(), 1080 )
  FCTEST_CHECK_EQUAL( mixed_triangles%size(), 1212 )
  FCTEST_CHECK( mixed_quads%size() > 0 )
  FCTEST_CHECK( mixed_triangles%size() > 0 )
  element_type = mixed_quads%element_type()
  FCTEST_CHECK_EQUAL( element_type%name(), "Quadrilateral" )
  call element_type%final()
  element_type = mixed_triangles%element_type()
  FCTEST_CHECK_EQUAL( element_type%name(), "Triangle" )
  call element_type%final()

  call mixed_quads%final()
  call mixed_triangles%final()
  call mixed_cells%final()
  call mixed_nodes%final()
  call mixed_mesh%final()
  call distribution%final()
  call partitioner%final()
  call grid%final()
END_TEST

TEST( test_mesh_grid_partitioner_config )
implicit none

  type(atlas_Grid) :: grid
  type(atlas_Partitioner) :: partitioner
  type(atlas_Config) :: config
  type(atlas_Mesh) :: triangular_mesh
  type(atlas_mesh_Nodes) :: triangular_nodes
  type(atlas_mesh_Cells) :: triangular_cells
  type(atlas_Elements) :: triangular_quads, triangular_triangles
  type(atlas_ElementType) :: element_type

  call atlas_log%info( "--- test_mesh_grid_partitioner_config" )

  grid = atlas_Grid("N16")
  partitioner = atlas_Partitioner("serial")
  config = atlas_Config()
  call config%set("triangulate", .true.)

  triangular_mesh = atlas_Mesh(grid, partitioner, config)

  triangular_nodes = triangular_mesh%nodes()
  FCTEST_CHECK( triangular_nodes%size() > 0 )

  triangular_cells = triangular_mesh%cells()
  FCTEST_CHECK_EQUAL( triangular_cells%nb_types(), 2 )
  triangular_quads = triangular_cells%elements(1)
  triangular_triangles = triangular_cells%elements(2)
  FCTEST_CHECK_EQUAL( triangular_nodes%size(), 1720 )
  FCTEST_CHECK_EQUAL( triangular_triangles%size(), 3372 )
  FCTEST_CHECK_EQUAL( triangular_quads%size(), 0 )
  FCTEST_CHECK( triangular_triangles%size() > 0 )
  element_type = triangular_triangles%element_type()
  FCTEST_CHECK_EQUAL( element_type%name(), "Triangle" )
  call element_type%final()

  call triangular_quads%final()
  call triangular_triangles%final()
  call triangular_cells%final()
  call triangular_nodes%final()
  call triangular_mesh%final()
  call config%final()
  call partitioner%final()
  call grid%final()
END_TEST

TEST( test_matching_partitioner )
implicit none

  type(atlas_Grid) :: grid
  type(atlas_Partitioner) :: partitioner, matching_partitioner
  type(atlas_Mesh) :: mesh
  type(atlas_FunctionSpace) :: functionspace

  call atlas_log%info( "--- test_matching_partitioner" )

  grid = atlas_Grid("N16")
  partitioner = atlas_Partitioner("serial")
  mesh = atlas_Mesh(grid, partitioner)

  matching_partitioner = atlas_MatchingPartitioner(mesh)
  call matching_partitioner%final()
  functionspace = atlas_functionspace_StructuredColumns(grid, partitioner)
  matching_partitioner = atlas_MatchingPartitioner(functionspace)

  call matching_partitioner%final()
  call functionspace%final()
  call mesh%final()
  call partitioner%final()
  call grid%final()
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

