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

TESTSUITE(fctest_ifs_setup)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  use atlas_module
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use atlas_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_fv )
use atlas_module
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

  write(0,*) "test_fv"

  ! Create a new Reduced Gaussian Grid based on a nloen array
  call atlas_log%info("Creating grid")
  grid = atlas_ReducedGaussianGrid( nloen(1:32) )
  FCTEST_CHECK_EQUAL( grid%owners(), 1 )

  ! Grid distribution: all points belong to partition 1
  allocate( part(grid%size()) )
  part(:) = 1
  griddistribution = atlas_GridDistribution(part, part0=1)
  FCTEST_CHECK_EQUAL( griddistribution%owners(), 1 )

  ! Generate mesh with given grid and distribution
  meshgenerator = atlas_MeshGenerator()
  FCTEST_CHECK_EQUAL( meshgenerator%owners(), 1 )

  write(0,*) " + mesh = meshgenerator%generate(grid,griddistribution)"
  write(0,*) "mesh%is_null()",mesh%is_null()
  mesh = meshgenerator%generate(grid,griddistribution)
  write(0,*) " + call griddistribution%final()"
  call griddistribution%final()

  ! Generate nodes function-space, with a given halo_size
  write(0,*) " + nodes_fs = atlas_functionspace_NodeColumns(mesh,halo_size)"
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
  FCTEST_CHECK_EQUAL( node_to_node%rows(), 0 )
  FCTEST_CHECK_EQUAL( node_to_node%name(),"node")

  node_to_node = nodes%connectivity("node")
  node_to_edge = nodes%connectivity("edge")

  call node_to_node%final()
  call mesh%final()
  call grid%final()
  call nodes_fs%final()

  write(0,*) "end test_fv"

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

